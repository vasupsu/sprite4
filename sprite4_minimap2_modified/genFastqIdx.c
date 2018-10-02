#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

typedef struct {
	size_t offset;
    size_t offset2;
	int fileNo;
//    uint32_t kmerFreqCount[1048576+2];
} fileIndex;

//long chunkSize=(long)25500000000;//reads
int readSize=230;
int chunkSize=50000000;
int numChunks=1;
int numPEFiles=0;//0 - single ended; 1 - interleaved; 2 - separate paired end fastq files
size_t *fileSize=NULL;
char **fileNames = NULL;

void getFileIndex (int chunkNum, size_t sizePerChunk, int numFiles, size_t *fileSize, fileIndex *result)
{
	int currentFileInd=0;
	size_t bytesRemaining = sizePerChunk * chunkNum;
	result->fileNo=currentFileInd;
	result->offset=bytesRemaining;
}
int file2No=1;
size_t file2Ofs=0;
#define READ_SIZE_MAX 100000000
size_t findOfsOtherEnd (int fileNo, size_t offset, char *readName, FILE *fp)
{
	int i=0;
//	printf ("Begin findOfsOtherEnd %s\n", readName);
	size_t delta = fileSize[1]/numChunks/2;
	if (delta > READ_SIZE_MAX) delta = READ_SIZE_MAX;
	char *fileBuffer = (char *)malloc(2*delta);
	assert (fileBuffer != NULL);
	char curReadName[1000];

        assert (fp != NULL);
	fseek (fp, offset, SEEK_SET);
	fgets (fileBuffer, 1000, fp);
	fileBuffer [strlen(fileBuffer)-1]='\0';
//	printf ("file2:%s\nfile1:%s\n", fileBuffer, readName);
	if ((strlen(fileBuffer)==strlen(readName)) && (strlen(fileBuffer)>2))
	{
		int equal=1;
		for (i=0; i<strlen (readName); i++)
		{
			if (fileBuffer[i]!=readName[i])
			{
				equal=0;
				break;
			}
			if ((fileBuffer[i]==' ') || (fileBuffer[i]=='/')) break;
		}
		if (equal) return offset;
	}
	
	for (i=0; i<strlen (readName); i++)
	{
		if ((readName[i] == ' ') || (readName[i]=='/')) 
		{
			readName[i] = '\0';
			break;
		}
	}
	
	size_t startOfs=0,len=0;
	if (offset >= delta)
		startOfs = offset-delta;
	else
		startOfs = 0;	
	len = fileSize[1] - startOfs;

	if ((offset + delta) < fileSize[fileNo+1])
		len = offset + delta - startOfs;
//	printf ("File1 ofs %lu seekOfs %lu, len %lu\n", offset, startOfs, len);
        fseek(fp, startOfs, SEEK_SET);
	fread (fileBuffer, 1, len, fp);
	i=0;
	
	while ((i<len) && ((fileBuffer[i] != '\n') || (fileBuffer[i+1] != '+')))
	{
		i++;
	}
	
	assert (i<len);
//	fileBuffer[i+100]='\0';
//	printf ("fileBuffer *%s*\n", &fileBuffer[i]);
	i++;
	while (fileBuffer[i] != '\n') i++;
	i++;
	while (fileBuffer[i] != '\n') i++;
        i++;
	fileBuffer[i+100]='\0';
//	printf ("fileBuffer *%s*\n", &fileBuffer[i]);

	if (fileBuffer[i] != '@')
	{
		while (fileBuffer[i] != '\n') i++;i++;
		assert (fileBuffer[i]=='+');
		while (fileBuffer[i] != '\n') i++;i++;
		while (fileBuffer[i] != '\n') i++;i++;
	}
	assert ((fileBuffer[i]=='@') && (i < len));
	int readNameInd=0;
	size_t j=i;
	while (fileBuffer[j] != '\n')
	{
		curReadName[readNameInd++] = fileBuffer[j];
		j++;
		if ((fileBuffer[j] == ' ') || (fileBuffer[j] == '/'))
			break;
	}
//	assert (curReadName[readNameInd-1] == '2');
	curReadName[readNameInd]='\0';
//	printf ("curReadName *%s*\n", curReadName);

	while (i < len)
	{
//		if (offset == 65077766728)
//			printf ("\tcurReadName *%s* readName *%s*\n", curReadName, readName);
		int comp=strcmp(curReadName, readName);
		if (comp == 0) break;
		while (fileBuffer[i] != '\n') {
			i++; 
		}
		i++;
		while (fileBuffer[i] != '\n') {
			i++; 
		}
		i++;
		while (fileBuffer[i] != '\n') {
			i++;
		}	
		i++;
		while (fileBuffer[i] != '\n') {
			i++;
		}	
		i++;
		readNameInd=0;
		j=i;
		while (fileBuffer[j] != '\n')
		{
			curReadName[readNameInd++] = fileBuffer[j];
			j++;
			if ((fileBuffer[j] == ' ') || (fileBuffer[j] == '/'))
				break;
		}
//		assert (curReadName[readNameInd-1] == '2');
		curReadName[readNameInd]='\0';
	}
//	printf ("\t2-%s Pos %lu\n", curReadName, startOfs+i);
//	printf ("findOfsOtherEnd end\n");
	free (fileBuffer);
	return (startOfs+i);
}


int main(int argc, char **argv) {

	int maxReadSize=300;
	//FILE *fp;
	struct stat buf;
	static struct option long_options[] =
	{
//		{"interleaved", required_argument,       0, 'p'},
		{"numchunks", required_argument,       0, 'c'},
		{0, 0, 0, 0}
	};
	int c=0;
	while (1)
	{
		int option_index = 0;
		c = getopt_long (argc, argv, "c:", long_options, &option_index);
		if (c==-1) break;
		switch (c)
		{
			case 'c':
				printf ("numchunks: %d\n", atoi(optarg));
				numChunks = atoi(optarg);
				break;
			default:
				printf ("unrecognized option -%c\n", c);
				return 1;
		}
	}
	int numPosArgs = argc-optind;
	if ((numPosArgs > 2) || (numPosArgs == 0) || ((numPosArgs > 1) && numPEFiles))
	{
		fprintf (stderr, "\n");
		fprintf (stderr, "Usage: %s [options] <in1.fq> [in2.fq]\n", argv[0]);
		fprintf (stderr, "Options:\n\n");
//		fprintf (stderr, "      -p            indicate in1.fq is an interleaved paired-end file. in2.fq argument not required\n");
		fprintf (stderr, "      -c INT        number of chunks [1]\n");
		return 0;
	}
	if (numPosArgs == 2)
		numPEFiles = 2;
//	printf ("numPEFiles %d\n", numPEFiles);
	//char readFile[500];

	int i;
	argv = &argv[optind];
	argc-=optind;
//	printf ("Argc %d\n", argc);

	fileNames = (char **)malloc(argc * sizeof(char *));
	fileSize = (size_t *)malloc(argc * sizeof(size_t));
	assert (fileNames != NULL);
	assert (fileSize != NULL);
	for (i=0; i<argc; i++)
	{
		fileNames[i]=(char *)malloc(500);
		assert (fileNames[i] != NULL);
        	sprintf (fileNames[i], "%s", argv[i]);
	        int status = stat (fileNames[i], &buf);
        	assert (status == 0);
		fileSize[i]=buf.st_size;
		printf ("file %s Size %lu\n", fileNames[i], buf.st_size);
//		if (i <= numPEFiles)
//			if ((i % 2) == 0) continue;
	}
	size_t sizePerChunk = fileSize[0]/numChunks;
//	printf ("size %lu sizePerChunk %lu\n", fileSize[0], sizePerChunk);
	size_t j;
	
	char *readChunk = (char *)malloc(5000);

	long writeOfs = 0;
	//FILE *fpOp = NULL;
	char idxFile[200];
	sprintf (idxFile, "%s.idx", fileNames[0]);
	//long startOfs=0;
	FILE *fpOut = fopen(idxFile, "w");
	assert (fpOut != NULL);
	//int eof=0;
	//size_t curOfs=2;
	fileIndex fI;
	fI.fileNo=0;
	fI.offset=0;

	int flag=0;	
	char curReadName[100];
	FILE *fp = fopen (fileNames[0], "r");
	assert (fp != NULL);
	FILE *fp2 = NULL;
	if (numPEFiles == 2)
	{
		fp2 = fopen (fileNames[1], "r");
		assert (fp2 != NULL);
	}
	for (i=0; i<numChunks; i++)
	{
		flag=0;
		getFileIndex (i, sizePerChunk, argc-1, fileSize, &fI);

		if (i==0)
		{
			fI.offset2=0;
			fwrite (&fI, sizeof(fileIndex), 1, fpOut);
			continue;
		}
		fseek(fp, fI.offset, SEEK_SET);
		//curOfs=startOfs;
		size_t nread =fread (readChunk, 1, 3000, fp);

//		printf ("nread %lu\n", nread);
//		assert (nread > maxReadSize);
		int k=0;
		for (j=0; j<nread; j++)
		{
			if ((readChunk[j]=='\n') && (readChunk[j+1]=='+'))
			{
//				if (readChunk[j+2]=='\n')
				{
					flag=1;
					j=j+2;
					while (readChunk[j] != '\n')
						j++;
					j++;
				}
			}
                        k=0;
			if (flag && (readChunk[j]=='\n'))
			{
				j++;k++;
				if (readChunk[j] != '@') flag = 0;
				else break;
			}
		}
//		if ((2*k+30) < maxReadSize) maxReadSize = 2*k+30;
		assert(flag==1);
		writeOfs = fI.offset+j;
//		printf ("writeOfs %lu\n", writeOfs);
		fI.offset = writeOfs;
		int readInd=0;
		while (readChunk[j]!='\n')
			curReadName[readInd++]=readChunk[j++];
//		assert (curReadName[readInd-1]=='1');
		curReadName[readInd]='\0';
//		printf ("FIle %d curReadName %s\n", fI.fileNo, curReadName);	
		if (fI.fileNo < numPEFiles)
		{	
			size_t ofs2=findOfsOtherEnd (fI.fileNo, fI.offset, curReadName, fp2);
			fI.offset2 = ofs2;
		}
		else
		{
			fI.offset2 = 0;
		}

		fwrite (&fI, sizeof(fileIndex), 1, fpOut);
//		printf ("P%d - %d %ld\n", rank, count, startOfs+i);
		printf ("Chunk %d: File %d, Offset %lu,%lu\n", i, fI.fileNo, fI.offset, fI.offset2);
	}
	fclose (fp);
	if (fp2) fclose (fp2);

	fI.fileNo = 0;
	fI.offset = fileSize[0];
	if (numPEFiles == 2)
		fI.offset2 = fileSize[1];
	else
		fI.offset2 = 0;

	fwrite (&fI, sizeof(fileIndex), 1, fpOut);
	printf ("Chunk %d: File %d, Offset %lu,%lu\n", i, fI.fileNo, fI.offset, fI.offset2);
	free (readChunk);
	for (i=0; i<argc; i++)
	{
		free (fileNames[i]);
	}
	free(fileNames);
	free (fileSize);
	fclose (fpOut);
	return 0;
}
