#include <stdlib.h>
#include <sys/stat.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "kthread.h"
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kvec.h"
#include "khash.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;
extern int numTasks, rank;
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

typedef struct mm_idx_bucket_s {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

void mm_idxopt_init(mm_idxopt_t *opt)
{
	memset(opt, 0, sizeof(mm_idxopt_t));
	opt->k = 15, opt->w = 10, opt->flag = 0;
	opt->bucket_bits = 14;
	opt->mini_batch_size = 50000000;
	opt->batch_size = 4000000000ULL;
}

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	int i;
	if (mi == 0) return;
	for (i = 0; i < 1<<mi->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	free(mi->B); free(mi->S); free(mi);
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

void mm_idx_stat(const mm_idx_t *mi)
{
	int i, n = 0, n1 = 0;
	uint64_t sum = 0, len = 0;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	for (i = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k)) {
				sum += kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
				if (kh_key(h, k)&1) ++n1;
			}
	}
	fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf\n",
			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), n, 100.0*n1/n, (double)sum / n, (double)len / sum);
}

int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	for (i = st1; i < en1; ++i)
		seq[i - st1] = mm_seq4_get(mi->S, i);
	return en - st;
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>8>>mi->b<<1, &absent);
			assert(absent && j - start_a == n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	size_t file_size;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x>>8&mask].a;
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		#pragma omp critical (read_lock)
		{
			s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq, UINT64_MAX); // read a mini-batch
		}
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		#pragma omp critical (comp_lock)
		{
			for (i = 0; i < s->n_seq; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				if (t->l_seq > 0)
					mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, p->mi->flag&MM_I_HPC, &s->a);
				else if (mm_verbose >= 2)
					fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
				free(t->seq); free(t->name);
			}
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		#pragma omp critical (comp_lock)
		{
			mm_idx_add(p->mi, s->a.n, s->a.a);
		}
		kfree(0, s->a.a); free(s);
	}
    return 0;
}

mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp)) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post(pl.mi, n_threads);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int flag, int n_threads) // a simpler interface; deprecated
{
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	fp = mm_bseq_open(fn, 0);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, 14, flag, 1<<18, n_threads, UINT64_MAX);
	mm_bseq_close(fp);
	return mi;
}

mm_idx_t *mm_idx_str(int w, int k, int is_hpc, int bucket_bits, int n, const char **seq, const char **name)
{
	uint64_t sum_len = 0;
	mm128_v a = {0,0,0};
	mm_idx_t *mi;
	int i, flag = 0;
	if (n <= 0) return 0;
	for (i = 0; i < n; ++i) // get the total length
		sum_len += strlen(seq[i]);
	if (is_hpc) flag |= MM_I_HPC;
	if (name == 0) flag |= MM_I_NO_NAME;
	if (bucket_bits < 0) bucket_bits = 14;
	mi = mm_idx_init(w, k, bucket_bits, flag);
	mi->n_seq = n;
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, n, sizeof(mm_idx_seq_t)); // ->seq is allocated from km
	mi->S = (uint32_t*)calloc((sum_len + 7) / 8, 4);
	for (i = 0, sum_len = 0; i < n; ++i) {
		const char *s = seq[i];
		mm_idx_seq_t *p = &mi->seq[i];
		uint32_t j;
		if (name && name[i]) {
			p->name = (char*)kmalloc(mi->km, strlen(name[i]) + 1);
			strcpy(p->name, name[i]);
		}
		p->offset = sum_len;
		p->len = strlen(s);
		for (j = 0; j < p->len; ++j) {
			int c = seq_nt4_table[(uint8_t)s[j]];
			uint64_t o = sum_len + j;
			mm_seq4_set(mi->S, o, c);
		}
		sum_len += p->len;
		if (p->len > 0) {
			a.n = 0;
			mm_sketch(0, s, p->len, w, k, i, is_hpc, &a);
			mm_idx_add(mi, a.n, a.a);
		}
	}
	free(a.a);
	mm_idx_post(mi, 1);
	return mi;
}

/*************
 * index I/O *
 *************/

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[5];
	int i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 5, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		l = strlen(mi->seq[i].name);
		fwrite(&l, 1, 1, fp);
		fwrite(mi->seq[i].name, 1, l, fp);
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

mm_idx_t *mm_idx_load(FILE *fp)
{
#define ITERS 5
	int i;
	char magic[4];
	uint32_t x[5];
	uint64_t sum_len = 0;
	mm_idx_t *mi;

	if (ftell(fp) != 0) return 0;
	fseek (fp, 0, SEEK_END);
	size_t idxFileSize = ftell (fp);
	printf ("Index File Size %lu\n", idxFileSize);
	size_t sizeToRead = idxFileSize/numTasks/ITERS;
	assert (sizeToRead > 0);
	fseek (fp, sizeToRead*rank, SEEK_SET);
	char *buf = (char *)malloc (idxFileSize+1);
	if (buf == NULL)
	{
		printf ("mm_idx_load: Cannot allocate memory\n");
		return 0;
	}
//	if (rank == 0)
	double stime = realtime();
	for (i=0; i<ITERS; i++)
	{
		fseek (fp, sizeToRead*numTasks*i+sizeToRead*rank, SEEK_SET);
		assert (fread (buf+sizeToRead*numTasks*i+sizeToRead*rank, 1, sizeToRead, fp) == sizeToRead);
		printf ("[%d] Iter %d File read %lu bytes time %lf sec\n", rank, i, sizeToRead, realtime()-stime);
		if (numTasks > 1)
		{
			MPI_Allgather (MPI_IN_PLACE, sizeToRead, MPI_BYTE, buf+sizeToRead*numTasks*i, sizeToRead, MPI_BYTE, MPI_COMM_WORLD);
			printf ("[%d] Allgather time %lf sec\n", rank, realtime()-stime);
		}
		MPI_Barrier (MPI_COMM_WORLD);
	}
	size_t ofs = sizeToRead*numTasks*ITERS;
	sizeToRead = idxFileSize - sizeToRead*numTasks*ITERS;
	
	if (sizeToRead > 0)
	{
		fseek (fp, ofs, SEEK_SET);
		assert (fread (buf+ofs, 1, sizeToRead, fp) == sizeToRead);
	}
	printf ("[%d] File read %lu bytes time %lf sec\n", rank, sizeToRead, realtime()-stime);
	MPI_Barrier (MPI_COMM_WORLD);
/*	if (rank != 0)
	{
		double stime = realtime();
		assert (fread (buf, 1, idxFileSize, fp) == idxFileSize);
		printf ("[%d] File read time %lf sec\n", rank, realtime()-stime);
	}
	MPI_Barrier (MPI_COMM_WORLD);*/
//	fseek (fp, 0, SEEK_SET);

	size_t bufOfs=0;
//	if (fread(magic, 1, 4, fp) != 4) return 0;
	if ((bufOfs+4) >= idxFileSize) return 0;
	memcpy (magic, &buf[bufOfs], 4);
	bufOfs += 4;
	
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
//	if (fread(x2, 4, 5, fp) != 5) return 0;
	if ((bufOfs+20) >= idxFileSize) return 0;
	memcpy (x, &buf[bufOfs], 20);
/*	if (! ((x[0] == x2[0]) && (x[1] == x2[1]) && (x[2] == x2[2]) && (x[3] == x2[3]) && (x[4]==x2[4])))
	{
		printf ("x[0-5] != x2[0-5] %u,%u,%u,%u,%u vs %u,%u,%u,%u,%u\n", x[0], x[1], x[2], x[3], x[4], x2[0], x2[1], x2[2], x2[3], x2[4]);
		return 0;
	}*/
	bufOfs += 20;

	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
//		fread(&l2, 1, 1, fp);
		l = *((uint8_t *)&buf[bufOfs]);
		bufOfs++;
//		assert (l==l2);
		s->name = (char*)kmalloc(mi->km, l + 1);
//		fread(s->name, 1, l, fp);
		strncpy (s->name, &buf[bufOfs], l);
		s->name[l] = 0;
		bufOfs += l;
//		fread(&s->len, 4, 1, fp);
		s->len = *((uint32_t *)&buf[bufOfs]);
		bufOfs+=4;

//		if (rank==0)
//			printf ("%d - %s %u\n", i, s->name, s->len);
		s->offset = sum_len;
		sum_len += s->len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
//		fread(&b->n, 4, 1, fp);
		b->n = *((int32_t *)&buf[bufOfs]);
		bufOfs+=4;
		b->p = (uint64_t*)malloc(b->n * 8);
//		fread(b->p, 8, b->n, fp);
		memcpy (b->p, &buf[bufOfs], b->n * 8);
		bufOfs += b->n * 8;

//		fread(&size2, 4, 1, fp);
		size = *((uint32_t *)&buf[bufOfs]);
//		assert (size == size2);
		bufOfs += 4;

		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
//			fread(x2, 8, 2, fp);
			memcpy (x, &buf[bufOfs], 2 * 8);
			bufOfs += 2 * 8;
//			assert ((x[0]==x2[0]) && (x[1] == x2[1]));

			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
//		fread(mi->S, 4, (sum_len + 7) / 8, fp);
		memcpy (mi->S, &buf[bufOfs], (sum_len + 7) / 8 *4);
		bufOfs += (sum_len + 7) / 8 *4;
	}
	free (buf);
	return mi;
}
#ifdef USE_MPIO
mm_idx_t *mm_idx_load_mpi(MPI_File fp)
{
	int i;
	char magic[4];
	uint32_t x[5];
	uint64_t sum_len = 0;
	mm_idx_t *mi;
	MPI_Status status;

//	if (fread(magic, 1, 4, fp) != 4) return 0;
	MPI_File_read_all (fp, magic, 4, MPI_CHAR, &status);
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
//	if (fread(x, 4, 5, fp) != 5) return 0;
	MPI_File_read_all (fp, x, 5, MPI_UNSIGNED, &status);
	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
//		fread(&l, 1, 1, fp);
		MPI_File_read_all (fp, &l, 1, MPI_BYTE, &status);
		s->name = (char*)kmalloc(mi->km, l + 1);
//		fread(s->name, 1, l, fp);
		MPI_File_read_all (fp, s->name, l, MPI_CHAR, &status);
		s->name[l] = 0;
//		fread(&s->len, 4, 1, fp);
		MPI_File_read_all (fp, &s->len, 1, MPI_UNSIGNED, &status);
		s->offset = sum_len;
		sum_len += s->len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
//		fread(&b->n, 4, 1, fp);
		MPI_File_read_all (fp, &b->n, 1, MPI_INT, &status);
		b->p = (uint64_t*)malloc(b->n * 8);
//		fread(b->p, 8, b->n, fp);
		MPI_File_read_all (fp, b->p, b->n, MPI_UNSIGNED_LONG, &status);
//		fread(&size, 4, 1, fp);
		MPI_File_read_all (fp, &size, 1, MPI_UNSIGNED, &status);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
//			fread(x, 8, 2, fp);
			MPI_File_read_all (fp, x, 2, MPI_UNSIGNED_LONG, &status);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
//		fread(mi->S, 4, (sum_len + 7) / 8, fp);
		MPI_File_read_all (fp, mi->S, (sum_len + 7) / 8, MPI_UNSIGNED, &status);
	}
	return mi;
}
#endif

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	off_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {
#ifdef USE_MPIIO
		 MPI_File_open (MPI_COMM_WORLD, fn, MPI_MODE_RDONLY, MPI_INFO_NULL, &(r->fp_idx));
#else
		r->fp.idx = fopen(fn, "rb");
#endif
		r->idx_size = is_idx;
	} 
	else 
	{
#ifdef USE_MPIIO
		r->fp_seq = mm_bseq_open(fn, 0);
#else
		r->fp.seq = mm_bseq_open(fn, 0);
#endif
	}
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
#ifdef USE_MPIIO
	if (r->is_idx) MPI_File_close (&(r->fp_idx));
	else mm_bseq_close(r->fp_seq);
#else
	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
#endif
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
	mm_idx_t *mi;
	if (r->is_idx) {
#ifdef USE_MPIIO
		mi = mm_idx_load_mpi(r->fp_idx);
#else
		mi = mm_idx_load(r->fp.idx);
#endif
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else
#ifdef USE_MPIIO
		mi = mm_idx_gen(r->fp_seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
#else
		mi = mm_idx_gen(r->fp.seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
#endif
	if (mi) {
		if (r->fp_out) mm_idx_dump(r->fp_out, mi);
		++r->n_parts;
	}
	return mi;
}

int mm_idx_reader_eof(const mm_idx_reader_t *r) // TODO: in extremely rare cases, mm_bseq_eof() might not work
{
#ifdef USE_MPIIO
	return 1;
#else
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
#endif
}