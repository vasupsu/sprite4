#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

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
int NUM_CONTIGS=0;

typedef struct
{
	char contigStr[100];
	long start;
	long end;
}interval;
int NUM_INTERVALS=0;
interval *iList=NULL;

void extract (char *prefix)
{
	int contigNo=0, pos=0;

	char *vcf_line = (char *)malloc (200000000);
	assert (vcf_line != NULL);
	char vcf_contigstr[200];
	int vcf_pos, vcf_qual;
	char bed_contigstr[200];
	int bed_pos;
	char fileName[300];

	sprintf (fileName, "%s_combined_extracted.vcf", prefix);
	FILE *fpOut = fopen (fileName, "w");
	assert (fpOut != NULL);

	sprintf (fileName, "%s_combined.vcf", prefix);
	FILE *fpVcf = fopen (fileName, "r");
	assert (fpVcf != NULL);
	fgets (vcf_line, 200000000, fpVcf);
	while (!feof (fpVcf) && (vcf_line[0]=='#'))
	{
		fputs (vcf_line, fpOut);
		fgets (vcf_line, 200000000, fpVcf);
	}
	if (!feof (fpVcf))
	{
		sscanf (vcf_line, "%s\t%d", vcf_contigstr, &vcf_pos);
	}
	int intervalNo=0;
	for (contigNo=0; contigNo<NUM_CONTIGS; contigNo++)
	{
		printf ("Contig %s\n", fileList[contigNo].fileName);
		for (pos=1; pos<=fileList[contigNo].seqLength; pos++)
		{
			int valid=0;
/*			if ((contigNo==0) && (pos==73))
			{
				printf ("curInterval %s %ld %ld VCF %s %ld\n", iList[intervalNo].contigStr, iList[intervalNo].start, iList[intervalNo].end, vcf_contigstr, vcf_pos);
			}*/
			if ((intervalNo<NUM_INTERVALS) && !strcmp (iList[intervalNo].contigStr, fileList[contigNo].fileName) && (iList[intervalNo].start <= pos) && (iList[intervalNo].end >= pos))
			{
				valid=1;
				if (iList[intervalNo].end == pos) 
				{
//					printf ("incrementing interval(%d:%ld-%ld) contigNo %d\n", intervalNo, iList[intervalNo].start, iList[intervalNo].end, contigNo);
					intervalNo++;
				}
			}
			if (!feof (fpVcf) && !strcmp (vcf_contigstr, fileList[contigNo].fileName) && (vcf_pos==pos))
			{
				if (valid)
					fputs (vcf_line, fpOut);
				while (!feof (fpVcf) && !strcmp (vcf_contigstr, fileList[contigNo].fileName) && (vcf_pos==pos))
				{
					fgets (vcf_line, 200000000, fpVcf);
					if (!feof (fpVcf))
						sscanf (vcf_line, "%s\t%d", vcf_contigstr, &vcf_pos);
				}
			}
		}
	}
	
	fclose (fpVcf);
	fclose (fpOut);

	free (vcf_line);
}

int main(int argc, char **argv)
{
	FILE *fpParsnip, *fpStrelka;
	if (argc != 4)
	{
		printf ("Usage: a.out <ref.fa.fai> <vcf prefix> <region.bed>\n");
		return 1;
	}
	
	char *fai_file = argv[1];
	printf ("faiFile %s\n", fai_file);
	char line[1000];
	long length=0,start=0;
	int contigNo=0;
	char contigName[100];
	FILE *fp_fai = fopen(fai_file, "r");
	assert (fp_fai != NULL);
	fgets(line, 200, fp_fai);
	while (!feof(fp_fai))
	{
		NUM_CONTIGS++;
		fgets(line, 200, fp_fai);
	}
	fclose (fp_fai);
	printf ("NUM_CONTIGS %d\n", NUM_CONTIGS);
	fileList = (fileInfo *)malloc(NUM_CONTIGS * sizeof(fileInfo));
	assert (fileList != NULL);

	fp_fai = fopen(fai_file, "r");
	assert (fp_fai != NULL);
	fgets(line, 200, fp_fai);
	while (!feof(fp_fai))
	{
		sscanf(line, "%s\t%ld\t%ld", contigName, &length, &start);
		fileList[contigNo].contigNo=contigNo;
		int fLen = strlen(contigName);
		fileList[contigNo].fileName = (char *)malloc(fLen+1);
		strcpy(fileList[contigNo].fileName, contigName);
		fileList[contigNo].seqLength = length;
		fileList[contigNo].REF_START = start;
		contigNo++;

		fgets(line, 200, fp_fai);
	}
	fclose (fp_fai);

//	printf ("%s
	FILE *fp_bed = fopen(argv[3], "r");
	assert (fp_bed != NULL);
	fgets(line, 200, fp_bed);
	while (!feof(fp_bed))
	{
		NUM_INTERVALS++;
		fgets(line, 200, fp_bed);
	}
	fclose (fp_bed);
	printf ("NUM_INTERVALS %d\n", NUM_INTERVALS);
	fp_bed = fopen(argv[3], "r");
	assert (fp_bed != NULL);
	fgets(line, 200, fp_bed);
	iList = (interval *)malloc (NUM_INTERVALS * sizeof(interval));
	assert (iList != NULL);
	int intervalNo=0;
	while (!feof(fp_bed))
	{
		sscanf (line, "%s\t%ld\t%ld\n", iList[intervalNo].contigStr, &iList[intervalNo].start, &iList[intervalNo].end);
		intervalNo++;
		fgets(line, 200, fp_bed);
	}
	
	fclose (fp_bed);
	extract(argv[2]);
	free (fileList);
	free (iList);
	return 0;
}
