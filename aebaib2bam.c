#include <stdio.h>
#include <ctype.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <getopt.h>
#include <stdint.h>
#include <sys/stat.h>
#include <assert.h>
#include "sam.h"
#include "hts.h"
#include "bgzf.h"
#include "kstring.h"
//#include "bam.h"
#define bam_reg2bin(b,e) hts_reg2bin((b),(e), 14, 5)
uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar) { return c->pos + (c->n_cigar? bam_cigar2rlen(c->n_cigar, cigar) : 1); }

char sam_prefix[500];
int NUM_CONTIGS=93;
char fai_file[200];
char bam_file[200];
char fileName[200];
char fasta_file[200];

int nextPos=-1;
int chunkSize = 262144*16;
int numThreads = 1;
FILE *fp_fasta;
int totalAlnFiles = 0;
int numtasks=1, rank=0;
int WINDOW_SIZE=50000, rSize=50000;
#define SEGMENT_SIZE 1000000
int maxBufferRecords=1000000;
typedef struct {
	int contigNo;
	int start;
	int end;
	int count;
}segmentInfo;
segmentInfo *si=NULL;
//FILTERS
/*
 * MQ - (avg_mapq >= 20)
 * BQ - (avg_snp_base_quality >= 10)
 * DP - (pos_table[ind].count[4] > 1)
 * SB - (r1>.1)
 * AAF - (maxPer > 20)
 * AAC - (max_per > 1)
 */ 
int MQ=20, DP=1, AAF=20, AAC=1, BQ=10;
float SB=.1;
int max_out_files=500;
int numOutChunks, numMaxChunks, maxChunkSize;

char **ref_chars=NULL;
typedef struct ref_pos
{
   uint32_t count[5];   
   uint32_t readPos[4];
   uint16_t fwdCount;//#alt alleles in fwd reads
   uint16_t fwdTotal;//#total fwd reads
   uint16_t revCount;
   uint16_t revTotal;
   uint32_t mapQ;
   uint16_t bqSnp;
}ref_pos;
typedef struct 
{
	int contigNo;
	char *fileName;
	long numRecords;
	long seqLength;
	long numAeb;
	long numAib;
        long REF_START;
}fileInfo;
fileInfo *fileList=NULL;

typedef struct
{
        uint32_t pos;
        uint8_t seq[60];
	uint8_t quals[120];
	uint16_t flag;
        uint8_t qual;
        uint8_t  matchLen;
}fullRec;

typedef struct
{
        uint32_t pos;
	uint16_t flag;
        uint8_t seq[60];
	uint8_t quals[120];
        uint8_t qual;
        uint8_t n_cigar;
        uint16_t cigar[10];
}otherRec;

typedef struct
{
 	fullRec *fRecs;
	otherRec *oRecs;
	int numFullRecs;
	int numOtherRecs;
}window;

long maxSeqLen=0;
int startContig=-1, endContig=-1;
int *procAssign=NULL;
int *startPos=NULL;
int *contigPosCount=NULL;

static int sizecmpfn(const void *a, const void *b) {
    fileInfo *av = (fileInfo *) a;
    fileInfo *bv = (fileInfo *) b;
    return av->numRecords-bv->numRecords;
}
void readFasta (int contigNo, FILE *fp_fasta, int len)
{
	long ind=0;
	fseek(fp_fasta, fileList[contigNo].REF_START, SEEK_SET);
	int lenToRead = len + (len/50)+100;
/*	if (contigStr==NULL)
	{
		printf ("Alloc fastaIp Str %ld\n", maxSeqLen + (maxSeqLen/50)+110);
	        contigStr = (char *)malloc(maxSeqLen + (maxSeqLen/50)+110);
		assert (contigStr != NULL);
	}*/
	char *contigStr = (char *)malloc(lenToRead+10);
        assert (contigStr != NULL);
        fread (contigStr, 1, lenToRead, fp_fasta);
	fclose (fp_fasta);
	int curIndex=0;
        for (ind=0; ind<len; ind++)
        {
	    char ch=contigStr[curIndex++];
            if (ch=='\n')
            {
                ind--;
                continue;
            }
            if (islower(ch))
                ch=toupper(ch);
            ref_chars[contigNo][ind]=ch;
        }
	free (contigStr);
}

int getSeqLen (otherRec *oR)
{
	int i=0, seqLen=0;
	for (i=0; i<oR->n_cigar; i++)
	{
		uint32_t len = oR->cigar[i]>>4;
		seqLen += len;
	}
	return seqLen;
}
void aeb2bam (bam1_t *b, fullRec *fR, int curTid, int segNo, int recNo)
{
	bam1_core_t *c = &b->core;
	uint8_t *curData = b->data;
	assert (curTid != -1);
	c->tid = curTid;
	c->pos = fR->pos-1;
	char rname[50];
	sprintf (rname, "R-%d-%d-%d", curTid, segNo, recNo);
	c->l_qname = strlen(rname)+1;
	uint64_t recordNo = curTid*numMaxChunks + segNo;
	recordNo = (recordNo << 32) | recNo;
	sprintf ((char *)curData, "%s", rname);
	curData[c->l_qname-1]=0; 
	
	c->n_cigar = 1;
	uint32_t cigar = fR->matchLen << 4;
	*((uint32_t *)(curData+c->l_qname)) = cigar;
	
	c->l_qseq = fR->matchLen;
	memcpy (&curData[c->l_qname+4], fR->seq, (c->l_qseq+1)>>1);
	memcpy (&curData[c->l_qname+4+((c->l_qseq+1)>>1)], fR->quals, c->l_qseq);
	
	
//	b->l_aux=0;//1?
	curData[c->l_qname+4+((c->l_qseq+1)>>1)+c->l_qseq]=0;
	
	b->l_data  = c->l_qname+4+((c->l_qseq+1)>>1)+c->l_qseq/*+1*/;
	b->m_data = b->l_data;
	b->m_data = kroundup32(b->m_data);
//	b->data = curData;
	
	
	uint32_t endPos = bam_calend(c, (uint32_t *)(curData+1));
	c->bin = bam_reg2bin(c->pos, endPos);
	c->qual = fR->qual;
	c->flag = fR->flag;
	c->mtid = -1;
	c->mpos = -1;
	c->isize = 0;
} 

void aib2bam (bam1_t *b, otherRec *oR, int curTid, int segNo, int recNo) 
{
	uint32_t cigarConv[6]={0,1,2,3,4,5};
	bam1_core_t *c = &b->core;
	uint8_t *curData = b->data;
	assert (curData != NULL);
	assert (curTid != -1);
	c->tid = curTid;
	c->pos = oR->pos-1;
	int i=0;
	
	char rname[50];
	sprintf (rname, "R-%d-%d-%d", curTid, segNo, recNo);
	c->l_qname = strlen(rname)+1;
	uint64_t recordNo = curTid*numMaxChunks + segNo;
	recordNo = (recordNo << 32) | recNo;
	sprintf ((char *)curData, "%s", rname);
	curData[c->l_qname-1]=0;
	
	c->n_cigar = oR->n_cigar;
	uint32_t *cigars = (uint32_t *)(curData+c->l_qname);
	c->l_qseq = 0;
	for (i=0; i<c->n_cigar; i++)
	{
	    assert ((oR->cigar[i] & 0xF) < 6);
	    cigars[i]=(oR->cigar[i] & 0xFFF0u) | cigarConv[oR->cigar[i] & 0xFu];
	    char c1="MIDNSH"[oR->cigar[i] & 0xF];
            assert (c1 != 'N');
	    uint32_t len = (oR->cigar[i])>>4;
	    if ((c1=='M') || (c1=='I') || (c1=='S')) c->l_qseq+=len;
  	    if ((segNo==3) && (recNo == 3))
	    {
		printf ("%s - %u%c,", rname, len, c1);
	    }
	}
  	    if ((segNo==3) && (recNo == 3))
	    	printf ("\n");
	memcpy (&curData[c->l_qname+(c->n_cigar*4)], oR->seq, (c->l_qseq+1)>>1);
	memcpy (&curData[c->l_qname+(c->n_cigar*4)+((c->l_qseq+1)>>1)], oR->quals, c->l_qseq);
	
	
//	b->l_aux=0;//1?
	curData[c->l_qname+(c->n_cigar*4)+((c->l_qseq+1)>>1)+c->l_qseq]=0;
	
	b->l_data = c->l_qname+(c->n_cigar*4)+((c->l_qseq+1)>>1)+c->l_qseq/*+1*/;
	b->m_data = b->l_data;
	b->m_data = kroundup32(b->m_data);
//	b->data = curData;
	
	
	uint32_t endPos = bam_calend(c, (uint32_t *)(curData+1));
	c->bin = bam_reg2bin(c->pos, endPos);
	c->qual = oR->qual;
	c->flag = oR->flag;
	c->mtid = -1;
	c->mpos = -1;
	c->isize = 0;
}

static int numComp (const void *p1, const void *p2)
{       
	int p = *(int *)p1;
	int q = *(int *)p2;
	if (p > q) return 1;
	else if (p < q) return -1;
	else return 0;
}

size_t writeBam (int segNo, bam_hdr_t *hdr, samFile* fpOut, fullRec *FR, otherRec *OR, int numAebs, int numAibs)
{
	int segStart = si[segNo].start;
	int segEnd = si[segNo].end;
	int i=0;
	size_t max_k = 0, k = 0;
	uint8_t *curData = (uint8_t *)calloc (10000, sizeof(uint8_t));
        assert (curData != NULL);
	bam1_t *b1=NULL;
        b1 = (bam1_t*)calloc(1, sizeof(bam1_t));
        assert (b1 != NULL);
        b1->data = curData;
	int curAeb=0, curAib=0;
	int nextOfs = rSize;
	char fName[300];
	while ((curAeb < numAebs) || (curAib < numAibs))
	{
		int curPos=0;
		fullRec *curFr = NULL;
		otherRec *curOr = NULL;
		if ((curAeb < numAebs) && ( ((curAib < numAibs) && (FR[curAeb].pos <= OR[curAib].pos)) || (curAib == numAibs) ))//processAEB rec
		{
			curPos = FR[curAeb].pos;
			curFr = &FR[curAeb];
			curAeb++;
		}
		else if ((curAib < numAibs) && ( ((curAeb < numAebs) && (OR[curAib].pos < FR[curAeb].pos)) || (curAeb == numAebs)) )//process AIB rec
		{
			curPos = OR[curAib].pos;
			curOr = &OR[curAib];
			curAib++;
		}
/*			if (k == max_k) {
				size_t old_max = max_k;
				max_k = max_k? max_k<<1 : 100;
//				printf ("Realloc  bytes (%lu*%lu)\n",  max_k, sizeof(bam1_t));
				buf = (bam1_t**)realloc(buf, max_k*sizeof(bam1_t *));
				assert (buf != NULL);
				memset(buf + old_max, 0, sizeof(bam1_t*) * (max_k - old_max));
			}*/
			if (curFr != NULL)
			{
				aeb2bam (b1, curFr, si[segNo].contigNo, si[segNo].start/maxChunkSize, curAeb+curAib);
			}
			else
			{
				aib2bam (b1, curOr, si[segNo].contigNo, si[segNo].start/maxChunkSize, curAeb+curAib);
			}
                        sam_write1(fpOut, hdr, b1);
                        k++;
//			printf ("curPos %d curIndelPos %d k %d\n", curPos, curIndelPos, k);
	}
	free (curData);
	return k;
}

void processSegment (int segNo, bam_hdr_t *hdr, samFile* fpOut, size_t *numBamRecs)
{
	struct stat sbuf;
//	printf ("[%d]process segment %d Contig %d:%d-%d Count %d\n", rank, segNo, si[segNo].contigNo, si[segNo].start, si[segNo].end, si[segNo].count);
	uint8_t *curData = (uint8_t *)calloc (10000, sizeof(uint8_t));
	assert (curData != NULL);

	int indelPosition=0, curIndelPos=0;

	int i=0;
	bam1_core_t *c;
	char *data;
	assert ((si[segNo].start % maxChunkSize) == 0);
	size_t max_k = 0, k = 0;
	bam1_t *b1=NULL;
	b1 = (bam1_t*)calloc(1, sizeof(bam1_t));
	assert (b1 != NULL);
	b1->data = curData;

	fullRec *FR = NULL;
	otherRec *OR = NULL;
	char fName[500];
	sprintf (fName, "%s/C%d_%d_sorted.aeb", sam_prefix, si[segNo].contigNo, si[segNo].start/maxChunkSize);
	int res = stat (fName, &sbuf);
	int numAebs = 0;
	if (res == 0)
	{
		numAebs = (int)(sbuf.st_size/sizeof(fullRec));
		FILE *fpAeb = fopen (fName, "r");
		assert (fpAeb != NULL);
		FR = (fullRec *)malloc (sbuf.st_size);
		assert (FR != NULL);
		fread (FR, sizeof(fullRec), numAebs, fpAeb);
		fclose (fpAeb);
	}
	sprintf (fName, "%s/C%d_%d_sorted.aib", sam_prefix, si[segNo].contigNo, si[segNo].start/maxChunkSize);
	res = stat (fName, &sbuf);
	int numAibs = 0;
	if (res == 0)
	{
		numAibs = (int)(sbuf.st_size/sizeof(otherRec));
		OR = (otherRec *)malloc (sbuf.st_size);
		assert (OR != NULL);
		FILE *fpAib = fopen (fName, "r");
		assert (fpAib != NULL);
		fread (OR, sizeof(otherRec), numAibs, fpAib);
		fclose (fpAib);
	}
//	printf ("segNo %d (%d,%d [%d,%d]) numAeb %d, numAib %d\n", segNo, si[segNo].contigNo, si[segNo].start/maxChunkSize, si[segNo].start, si[segNo].end, numAebs, numAibs);
	int curAeb=0, curAib=0;
	int *posList = NULL, nPos = 0, mPos = 0;
	if ((numAebs > 0) || (numAibs > 0))
	{
		*numBamRecs = writeBam (segNo, hdr, fpOut, FR, OR, numAebs, numAibs);
		printf ("Seg %d nBAMRecs %lu\n", segNo, *numBamRecs);
		*numBamRecs = nPos;
	}
	else
		*numBamRecs = 0;
/*	for (k = 0; k < max_k; ++k) {
		if (!buf[k]) continue;
		free(buf[k]->data);
		free(buf[k]);
	}*/
	if (FR)
		free (FR);
	if (OR)
		free (OR);
	free (curData);
//	free(buf);
/*	hts_idx_destroy (idx);
	hts_close (fp);
	hts_itr_destroy (iter);
	bam_hdr_destroy(hdr);
	bam_destroy1(b);*/
	
}

void split_fasta (int numChunks, fileInfo *anns, int numContigs, int *numOutChunks, int *numMaxChunks, int *maxChunkSize)
{
        int32_t maxContigSize = 0;
        int i=0;
        for (i=0; i<numContigs; i++)
        {
                if (anns[i].seqLength > maxContigSize) maxContigSize = anns[i].seqLength;
        }
        int curChunks = 0, curContigSize=maxContigSize;
        while (curChunks < numChunks)
        {
                curChunks = 0;
                for (i=0; i<numContigs; i++)
                {
                        curChunks += (anns[i].seqLength/curContigSize)+(anns[i].seqLength%curContigSize > 0);
                }
                if (rank == 0)
                        printf ("curChunkSize %d --> numChunks %d\n", curContigSize, curChunks);
                if (curChunks < numChunks)
                        curContigSize = curContigSize / 2;
        }
        if (curContigSize < maxContigSize)
                curContigSize = curContigSize * 2;
        int maxChunksPerContig = 0;
        curChunks = 0;
        for (i=0; i<numContigs; i++)
        {
                int curChunksPerContig = (anns[i].seqLength/curContigSize)+(anns[i].seqLength%curContigSize > 0);
                curChunks += curChunksPerContig;
                if (curChunksPerContig > maxChunksPerContig) maxChunksPerContig = curChunksPerContig;
        }
        *numOutChunks = curChunks;
        *numMaxChunks = maxChunksPerContig;
        *maxChunkSize = curContigSize;
}
void copyBamFile (samFile* fpOut, char *infile, int *first, int tid)
{       
        BGZF *fp=bgzf_open(infile, "r");
        assert (fp!=NULL);
        bam_hdr_t *hdr = bam_hdr_read(fp);
        bam1_t *b = (bam1_t *)calloc(1, sizeof(bam1_t));
        printf ("[%d]\tWriting %s\n", tid, infile);
        if (*first)
        {       
                *first=0;
                sam_hdr_write(fpOut, hdr);
        }
        while (bam_read1(fp, b) >= 0)
        {       
                sam_write1(fpOut, hdr, b);
        }       
        bgzf_close (fp);
        bam_hdr_destroy(hdr);
        bam_destroy1(b);
//	int stat=unlink(infile);
//	if (stat==-1)
//		printf ("Cannot remove BAM file %s. Error code %d\n", infile, errno);
}
int main(int argc, char **argv)
{
        struct stat buf;
        int i=0,j=0;
	static struct option long_options[] =
	{
		{"OUT", required_argument,       0, 'o'},
		{0, 0, 0, 0}
	};
	int c=0;
	while (1)
        {
		int option_index = 0;
		c = getopt_long (argc, argv, "o:", long_options, &option_index);
		if (c==-1) break;
		switch (c)
		{
			case 'o':
				max_out_files=(int)atoi(optarg);
				break;
			default:
				printf ("unrecognized option -%c\n", c);
				return 1;
		}
	}
        if (optind + 1 >= argc || optind + 3 != argc) 
        {
		fprintf (stderr, "\n");
                fprintf (stderr, "Usage: aebaib2bam [options] <path to fasta file> <Sorted AEB, AIB prefix> <output BAM filename>\n");
                fprintf (stderr, "Options:\n\n");
		fprintf (stderr, "	--OUT=INT     maximum number of FASTA regions [500]\n");
                return 0;
        }
        strcpy(fasta_file, argv[optind]);
        sprintf (fai_file, "%s.fai", fasta_file);
	FILE *fp_fai = fopen(fai_file, "r");
        strcpy (sam_prefix, argv[optind+1]);
        strcpy (bam_file, argv[optind+2]);

        char line[200];
        long length=0,start=0;

        NUM_CONTIGS=0;
        fgets(line, 200, fp_fai);
        while (!feof(fp_fai))
        {
                NUM_CONTIGS++;
                fgets(line, 200, fp_fai);
        }
        fclose(fp_fai);
        fileList = (fileInfo *)malloc(NUM_CONTIGS * sizeof(fileInfo));
        fp_fai = fopen(fai_file, "r");
        fgets(line, 200, fp_fai);
        int contigNo=0;
        char contigName[100];
        char fileName[500];
        long totalRecords=0;

	size_t genomeSize=0;
        while (!feof(fp_fai))
        {
                sscanf(line, "%s\t%ld\t%ld", contigName, &length, &start);
		fileList[contigNo].contigNo=contigNo;
                int fLen = strlen(contigName);
                fileList[contigNo].fileName = (char *)malloc(fLen+1);
                strcpy(fileList[contigNo].fileName, contigName);
                fileList[contigNo].numRecords = 0;
                fileList[contigNo].seqLength = length;
		genomeSize += length;
                fileList[contigNo].REF_START = start;
                if (maxSeqLen < length)
                        maxSeqLen = length;
                contigNo++;
                fgets(line, 200, fp_fai);
        }
        fclose(fp_fai);
	split_fasta (max_out_files, fileList, NUM_CONTIGS, &numOutChunks, &numMaxChunks, &maxChunkSize);
	if (rank == 0)
                printf ("NUM_CONTIGS %d numOutChunks %d numMaxChunks %d maxChunkSize %d\n", NUM_CONTIGS, numOutChunks, numMaxChunks, maxChunkSize);
	printf ("[%d]reading FASTA.. \n", rank);
	
        int totalSegments=0;
        int *numSegmentsPerContig = (int *)calloc (NUM_CONTIGS, sizeof(int));
        assert (numSegmentsPerContig != NULL);

	for (i=0; i<NUM_CONTIGS; i++)
	{
                numSegmentsPerContig[i] = fileList[i].seqLength/maxChunkSize;
                if ((fileList[i].seqLength % maxChunkSize) != 0) numSegmentsPerContig[i]++;
                totalSegments += numSegmentsPerContig[i];
	}
	printf ("Done \n totalSegments %d\n", totalSegments);
	si = (segmentInfo *)malloc (totalSegments * sizeof (segmentInfo));
        int curSegNo = 0;
        for (i=0; i<NUM_CONTIGS; i++)
        {
                for (j=0; j<numSegmentsPerContig[i]; j++)
                {
                        si[curSegNo].contigNo = i;
                        si[curSegNo].start = j*maxChunkSize;
                        if (j==(numSegmentsPerContig[i]-1))
                                si[curSegNo].end = fileList[i].seqLength;
                        else
                                si[curSegNo].end = (j+1)*maxChunkSize;
                        si[curSegNo].count = 0;
//                        if (i==1)
//                                printf ("segNo %d start %d end %d\n", curSegNo, si[curSegNo].start, si[curSegNo].end);
                        curSegNo++;
                }
        }
	free (numSegmentsPerContig);

	int tid=0;
      	struct timeval st;
        gettimeofday(&st, NULL);
	kstring_t str;
	str.l = str.m = 0; str.s = 0;
	char tStr[200];
	for (i=0; i<NUM_CONTIGS; i++)
	{
		sprintf (tStr, "@SQ\tSN:%s\tLN:%d\n", fileList[i].fileName, fileList[i].seqLength);
		kputsn (tStr, strlen(tStr), &str);
	}
//	sprintf (tStr, "@PG\tID:parsnip VN:4 CL:parsnip_aebaib %s %s %s\n", fasta_file, sam_prefix, vcf_file);
//	kputsn (tStr, strlen(tStr), &str);
//	printf ("*%s*%d\n", str.s, str.l);
	samFile* fpOut = sam_open(bam_file, "wb1");
	assert (fpOut != NULL);
	bam_hdr_t *hdr = sam_hdr_parse(str.l, str.s);
	sam_hdr_write(fpOut, hdr);
	hts_set_threads(fpOut, 1);
//	printf ("hdr->n_targets %d, name %s\n", hdr->n_targets, hdr->target_name[0]);
	size_t *sizes = (size_t *)calloc ((totalSegments + 2), sizeof(size_t));
	assert (sizes != NULL);///////////////////vasu
	for (i=0; i<totalSegments; i++)
	{
		processSegment (i, hdr, fpOut, &sizes[i+2]);
	}

        sam_close(fpOut);
	free (si);
	return 0;
}
