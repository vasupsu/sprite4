#include <stdio.h>
#include <ctype.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OMP
#include <omp.h>
#endif
#include <math.h>
#include <getopt.h>
#include <stdint.h>
#include <sys/stat.h>
#include <assert.h>
#include "sam.h"
#include "hts.h"
#include "bgzf.h"
#include "kstring.h"

#define bam_reg2bin(b,e) hts_reg2bin((b),(e), 14, 5)
uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar) { return c->pos + (c->n_cigar? bam_cigar2rlen(c->n_cigar, cigar) : 1); }
unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void aeb2bam (bam1_t *b, int pos, int matchLen, char *charseq, int curTid, bam_hdr_t *hdr, samFile *fpOut)
{
	uint8_t seq[60], quals[500];
	uint8_t bToInt[5] = {1, 2, 4, 8, 15};
        bam1_core_t *c = &b->core;
        uint8_t *curData = b->data;
        assert (curTid != -1);
        c->tid = curTid;
        c->pos = pos-1;
        char rname[50];
        sprintf (rname, "R-%d", curTid);
        c->l_qname = strlen(rname)+1;
        sprintf ((char *)curData, "%s", rname);
        curData[c->l_qname-1]=0;

        c->n_cigar = 1;
        uint32_t cigar = matchLen << 4;
        *((uint32_t *)(curData+c->l_qname)) = cigar;

        c->l_qseq = matchLen;
	int i=0;
	for (i=0; i<matchLen; i++)
	{
		int seqInd = i/2;
		int subInd = i&1;
		if (subInd==0)
			seq[seqInd] = bToInt[seq_nt4_table[(int)charseq[i]]]<<4;
		else
			seq[seqInd] |= bToInt[seq_nt4_table[(int)charseq[i]]];
		quals[i]=30;
	}
	
        memcpy (&curData[c->l_qname+4], seq, (c->l_qseq+1)>>1);
        memcpy (&curData[c->l_qname+4+((c->l_qseq+1)>>1)], quals, c->l_qseq);

        curData[c->l_qname+4+((c->l_qseq+1)>>1)+c->l_qseq]=0;

        b->l_data  = c->l_qname+4+((c->l_qseq+1)>>1)+c->l_qseq/*+1*/;
        b->m_data = b->l_data;
        b->m_data = kroundup32(b->m_data);

        uint32_t endPos = bam_calend(c, (uint32_t *)(curData+1));
        c->bin = bam_reg2bin(c->pos, endPos);
        c->qual = 60;
        c->flag = 0;
        c->mtid = -1;
        c->mpos = -1;
        c->isize = 0;

	sam_write1(fpOut, hdr, b);
}

void getReads (char *fasta_file, long *offsets, long *lengths, int numContigs, bam_hdr_t *hdr, samFile *fpOut)
{
	int i=0;
	char readStr[200];
	int readLen = 100;
	bam1_t *b1=NULL;
	b1 = (bam1_t*)calloc(1, sizeof(bam1_t));
	assert (b1 != NULL);
	uint8_t *curData = (uint8_t *)calloc (10000, sizeof(uint8_t));
	assert (curData != NULL);
	b1->data = curData;

	FILE *fp_fasta = fopen (fasta_file, "r");
	assert (fp_fasta != NULL);
	for (i=0; i<numContigs; i++)
	{
		fseek (fp_fasta, offsets[i], SEEK_SET);
		int start=0, startOfs=0, readOfs=0;
		char ch = fgetc(fp_fasta);
		if (ch != 'N') start=1;
		while (!start)
		{
			if (ch != '\n') startOfs++;
			ch = fgetc(fp_fasta);
			if ((ch != 'N') && (ch != '\n')) start=1;
		}
		readStr[0] = toupper(ch);
		readOfs=1;
		while (readOfs != readLen)
		{
			readStr[readOfs] = toupper(fgetc(fp_fasta));
			if (readStr[readOfs] != '\n') readOfs++;
		}
		readStr[readLen] = '\0';
//		printf ("Contig %d: %ld:%s\n", i, startOfs+1, readStr);
		aeb2bam (b1, (int)startOfs+1, readLen, readStr, i, hdr, fpOut);
	}
	fclose (fp_fasta);
	free (b1);
	free (curData);
}
int main(int argc, char **argv)
{
	if (argc != 3)
	{
		printf ("Usage: a.out <ref.fa> <out.bam>\n");
		return 1;
	}
	char line[1000];
	char fai_file[200], *bam_file, *fasta_file;
	fasta_file = argv[1];
	sprintf (fai_file, "%s.fai", fasta_file);
	bam_file = argv[2];
	
	FILE *fp_fai = fopen(fai_file, "r");
	assert (fp_fai != NULL);
	int NUM_CONTIGS=0;
	fgets(line, 200, fp_fai);
	while (!feof(fp_fai))
	{
		NUM_CONTIGS++;
		fgets(line, 200, fp_fai);
	}
	fclose(fp_fai);
//	printf ("# contigs: %d\n", NUM_CONTIGS);
	
	long *offsets = (long *)malloc (NUM_CONTIGS * sizeof(long));
	assert (offsets != NULL);
	long *lengths = (long *)malloc (NUM_CONTIGS * sizeof(long));
	assert (lengths != NULL);
	int contigNo=0;
	char contigName[100];
	long length=0,start=0;
	kstring_t str;
	char tStr[200];
	fp_fai = fopen(fai_file, "r");
	fgets(line, 200, fp_fai);
	while (!feof(fp_fai))
	{
		sscanf(line, "%s\t%ld\t%ld", contigName, &length, &start);
//		printf ("%s\t%ld\t%ld\n", contigName, length, start);
		offsets[contigNo] = start;
		lengths[contigNo] = length;
		sprintf (tStr, "@SQ\tSN:%s\tLN:%d\n", contigName, length);
		kputsn (tStr, strlen(tStr), &str);
		contigNo++;
		fgets(line, 200, fp_fai);
	}
	fclose(fp_fai);

	bam_hdr_t *hdr = sam_hdr_parse(str.l, str.s);
	samFile*  fpOut = sam_open(bam_file, "wb1");
	hts_set_threads(fpOut, 1);
	assert (fpOut != NULL);
	sam_hdr_write(fpOut, hdr);
	getReads (fasta_file, offsets, lengths, NUM_CONTIGS, hdr, fpOut);
	sam_close(fpOut);

	free (offsets);
	free (lengths);
	return 0;
}
