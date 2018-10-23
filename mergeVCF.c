#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

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

void combineStrelkaVcfs (char *prefix, int numFiles)
{
	char fileName[300];
	char *strelka_line = (char *)malloc (200000000);
	assert (strelka_line != NULL);
	sprintf (fileName, "%s_strelka2.vcf", prefix);
	int first=1, i=0;
	FILE *fpOut = fopen (fileName, "w");
	assert (fpOut != NULL);
	for (i=0; i<numFiles; i++)
	{
		sprintf (fileName, "%s.variants%d.vcf", prefix, i);
		FILE *fpIn = fopen (fileName, "r");
		if (fpIn == NULL) continue;
		fgets (strelka_line, 200000000, fpIn);
		while (!feof (fpIn))
		{
			if (strelka_line[0]=='#')
			{
				if (first)
					fputs (strelka_line, fpOut);
			}
			else
				fputs (strelka_line, fpOut);
			fgets (strelka_line, 200000000, fpIn);
		}
		first=0;
		fclose (fpIn);
	}
	fclose (fpOut);
	free (strelka_line);
}

void mergeParsnipStrelka (char *parsnip_prefix, char *varcall_prefix, char *out_prefix)
{
	int contigNo=0, pos=0;

	char *parsnip_line = (char *)malloc (2000);
	assert (parsnip_line != NULL);
	char *strelka_line = (char *)malloc (200000000);
	assert (strelka_line != NULL);
	char strelka_contigstr[200], posName[100], refBase[100], altBase[100], filter[200];
	int strelka_pos, strelka_qual;
	char parsnip_contigstr[200];
	int parsnip_pos;
	char fileName[300];

	sprintf (fileName, "%s", parsnip_prefix);
	FILE *fpParsnip = fopen (fileName, "r");
	assert (fpParsnip != NULL);
	fgets (parsnip_line, 2000, fpParsnip);
	while (!feof (fpParsnip) && (parsnip_line[0]=='#'))
	{
		fgets (parsnip_line, 2000, fpParsnip);
	}
	if (!feof (fpParsnip))
	{
		sscanf (parsnip_line, "%s\t%d", parsnip_contigstr, &parsnip_pos);
	}

	sprintf (fileName, "%s", out_prefix);
	FILE *fpOut = fopen (fileName, "w");
	assert (fpOut != NULL);

	sprintf (fileName, "%s", varcall_prefix);
	FILE *fpStrelka = fopen (fileName, "r");
	assert (fpStrelka != NULL);
	fgets (strelka_line, 200000000, fpStrelka);
	while (!feof (fpStrelka) && (strelka_line[0]=='#'))
	{
		fputs (strelka_line, fpOut);
		fgets (strelka_line, 200000000, fpStrelka);
	}
	if (!feof (fpStrelka))
	{
		sscanf (strelka_line, "%s\t%d\t%s\t%s\t%s\t%d\t%s", strelka_contigstr, &strelka_pos, posName, refBase, altBase, &strelka_qual, filter);
	}
	
	for (contigNo=0; contigNo<NUM_CONTIGS; contigNo++)
	{
		printf ("Contig %s\n", fileList[contigNo].fileName);
		for (pos=1; pos<=fileList[contigNo].seqLength; pos++)
		{
			int done=0;
			if (!feof (fpStrelka) && !strcmp (strelka_contigstr, fileList[contigNo].fileName) && (strelka_pos==pos))
			{
				if ((!strcmp (filter, "PASS")) && (done==0))
				{
					fputs (strelka_line, fpOut);
					done=1;
					fgets (strelka_line, 200000000, fpStrelka);
					if (!feof (fpStrelka))
						sscanf (strelka_line, "%s\t%d\t%s\t%s\t%s\t%d\t%s", strelka_contigstr, &strelka_pos, posName, refBase, altBase, &strelka_qual, filter);
				}
				while (!feof (fpStrelka) && !strcmp (strelka_contigstr, fileList[contigNo].fileName) && (strelka_pos==pos))
				{
					if (!strcmp (filter, "PASS") && (done==0))
					{
						done=1;
						fputs (strelka_line, fpOut);
					}
					fgets (strelka_line, 200000000, fpStrelka);
					if (!feof (fpStrelka))
						sscanf (strelka_line, "%s\t%d\t%s\t%s\t%s\t%d\t%s", strelka_contigstr, &strelka_pos, posName, refBase, altBase, &strelka_qual, filter);
				}
			}
			if (!feof (fpParsnip) && !strcmp (parsnip_contigstr, fileList[contigNo].fileName) && (parsnip_pos==pos))
			{
				if (!done)
				{
					fputs (parsnip_line, fpOut);
					done=1;
				}
				fgets (parsnip_line, 2000, fpParsnip);
				if (!feof (fpParsnip))
					sscanf (parsnip_line, "%s\t%d", parsnip_contigstr, &parsnip_pos);
			}
		}
	}
	
	fclose (fpStrelka);
	fclose (fpParsnip);
	fclose (fpOut);

	free (parsnip_line);
	free (strelka_line);
}

int main(int argc, char **argv)
{
	FILE *fpParsnip, *fpStrelka;
	if (argc != 5)
	{
		printf ("Usage: a.out <ref.fa.fai> <parsnip vcf> <varcall vcf> <out vcf>\n");
		return 1;
	}
//	int numFiles = atoi(argv[3]);
//	combineStrelkaVcfs (argv[2], numFiles);
	
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

	mergeParsnipStrelka(argv[2], argv[3], argv[4]);
	free (fileList);
	return 0;
}
