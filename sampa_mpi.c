#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OMP
#include <omp.h>
#endif
#include <stdint.h>
#include <sys/stat.h>
#include <assert.h>
#include <sys/time.h>

int bwa_P;
char sam_prefix[500];
int NUM_CONTIGS=93;
int NUM_FILES=0;
int max_out_files = 0;
int numOutChunks, numMaxChunks, maxChunkSize;
char fai_file[1000];
int numtasks=1, rank=0, numThreads = 1, rSize=50000;
typedef struct 
{
	char *fileName;
	long numRecords;
	long seqLength;
}fileInfo;
typedef struct
{       
        uint32_t pos;
        uint8_t seq[150];
        uint8_t quals[300];
        uint16_t flag;
        uint8_t qual;
        uint8_t  matchLen;
}fullRec;

typedef struct
{
        uint32_t pos;
        uint16_t flag;
        uint8_t seq[150];
        uint8_t quals[300];
        uint8_t qual;
        uint8_t n_cigar;
        uint16_t cigar[10];
}otherRec;

fileInfo *fileList=NULL;
long maxSeqLen=0;
int startContig=-1, endContig=-1;

static int sizecmpfn_aeb(const void *a, const void *b) {
    fullRec *av = (fullRec *) a;
    fullRec *bv = (fullRec *) b;
    return (av->pos - bv->pos);
}
static int sizecmpfn_aib(const void *a, const void *b) {
    otherRec *av = (otherRec *) a;
    otherRec *bv = (otherRec *) b;
    return (av->pos - bv->pos);
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
//		if (rank == 0)
//	                printf ("curChunkSize %d --> numChunks %d\n", curContigSize, curChunks);
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
long mpi_gatherv(void *lFr, size_t sendCount2, void *Fr2, size_t *recvCounts2, size_t *displ2, int root, int numtasks, int rank)
{
	int MAX_SIZE=1000000000, i=0;
	size_t j=0;
	struct timeval stime, etime;
	gettimeofday (&stime, NULL);
	for (i=0; i<numtasks; i++)
	{
		for (j=0; j<recvCounts2[i]; j+=MAX_SIZE)
		{
			int size_to_send = (j+MAX_SIZE > recvCounts2[i]) ? recvCounts2[i]-j : MAX_SIZE;
			if (rank == i) {
				if (rank == root) {
					memcpy ((void *)((char*)Fr2+displ2[i]+j), (void *)((char *)lFr+j), size_to_send);
				}
				else {
//					printf ("%d sending %dB to %d tag %d\n", rank, size_to_send, root, j);
					MPI_Send ((void *)((char *)lFr+j),size_to_send, MPI_UNSIGNED_CHAR, root, j, MPI_COMM_WORLD);
//					printf ("%d sent %dB to %d tag %d\n", rank, size_to_send, root, j);
				}
			}
			if (rank==root && rank!=i) {
//				printf ("%d recving %dB from %d tag %d\n", rank, size_to_send, i, j);
				MPI_Recv ((void *)((char*)Fr2+displ2[i]+j), size_to_send, MPI_UNSIGNED_CHAR, i, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//				printf ("%d recvd %dB from %d tag %d\n", rank, size_to_send, i, j);
			}
		}
	}
	gettimeofday (&etime, NULL);
	long elapsed = (etime.tv_sec * 1000000 + etime.tv_usec) - (stime.tv_sec * 1000000 + stime.tv_usec);
	return elapsed;
}
void write_histfile(char *sam_prefix, int curContigSegment, int contigNo, int segId, size_t *totSizes_aeb, size_t *totSizes_aib, fullRec *Fr2, otherRec *Or2)
{
	char fileName[500];
	sprintf (fileName, "%s0/C%d_%d_sorted.hist", sam_prefix, contigNo, segId);
	FILE *fp_hist = fopen (fileName, "wb");
	assert (fp_hist != NULL);
	int aebIdx = 0, aibIdx = 0;
	unsigned int minOfs = -1;
	if (aebIdx<totSizes_aeb[curContigSegment] && Fr2[aebIdx].pos<minOfs) {
		minOfs = Fr2[aebIdx].pos;
	}
	if (aibIdx<totSizes_aib[curContigSegment] && Or2[aibIdx].pos<minOfs) {
		minOfs = Or2[aibIdx].pos;
	}
	if (aebIdx<totSizes_aeb[curContigSegment] && minOfs==Fr2[aebIdx].pos)
		aebIdx++;
	else
		aibIdx++;
	fwrite(&minOfs, sizeof(int), 1, fp_hist);
	int curCount = 0;
	while (aebIdx<totSizes_aeb[curContigSegment] || aibIdx<totSizes_aib[curContigSegment]) //wip
	{
		for (curCount=0; curCount<rSize; curCount++) {
			if (aibIdx >= totSizes_aib[curContigSegment])
			{
				minOfs = Fr2[aebIdx].pos;
				aebIdx++;
			}
			else if (aebIdx >= totSizes_aeb[curContigSegment])
			{
				minOfs = Or2[aibIdx].pos;
				aibIdx++;
			}
			else
			{
				if (Fr2[aebIdx].pos <  Or2[aibIdx].pos)
				{
					minOfs = Fr2[aebIdx].pos;
					aebIdx++;
				}
				else
				{
					minOfs = Or2[aibIdx].pos;
					aibIdx++;
				}
			}
			if (aebIdx>=totSizes_aeb[curContigSegment] && aibIdx >= totSizes_aib[curContigSegment])
				break;
		}
		fwrite(&minOfs, sizeof(int), 1, fp_hist);
	}
	fclose (fp_hist);
}
int main(int argc, char **argv)
{
	struct stat buf;
	int i=0,j=0,k=0;
	static struct option long_options[] =
	{
		{"nthreads", required_argument,       0, 't'},
		{0, 0, 0, 0}
	};
	int c=0;
	while (1)
	{
		int option_index = 0;
		c = getopt_long (argc, argv, "t:", long_options, &option_index);
		if (c==-1) break;
		switch (c)
		{
			case 't':
				numThreads = atoi(optarg);
				break;
			default:
				printf ("unrecognized option -%c\n", c);
				return 1;
		}
	}
	if (optind + 4 != argc)
	{
		fprintf (stderr, "\n");
                fprintf (stderr, "Usage: sampa [options] <maxOutFiles> <ref.fasta.fai> <ain,aeb file prefix> <#MAP MPI processes>\n");
                fprintf (stderr, "Positional arguments:\n\n");
                fprintf (stderr, "maxOutFiles           Max number of output files per process\n");
                fprintf (stderr, "ref.fasta.fai         Absolute path of the fai file corresponding to the input fasta file\n");
                fprintf (stderr, "aeb,aib file prefix   prefix of the filename of aeb and aib files created by MAP step\n");
                fprintf (stderr, "#MAP MPI processes    Number of MPI tasks for Minimap2 run\n");
                fprintf (stderr, "\nOptions:\n\n");
                fprintf (stderr, "      -t INT          number of threads [1]\n");
		return 0;
	}
	max_out_files = atoi (argv[optind]);
	strcpy(fai_file, argv[optind+1]);
	FILE *fp_fai = fopen(fai_file, "r");
	assert(fp_fai != NULL);
	strcpy (sam_prefix, argv[optind+2]);
	bwa_P = atoi(argv[optind+3]);
	printf ("FAI: %s\nMAX_SEGS: %d\nPrefix path: %s\nMinimap2 tasks: %d\n", fai_file, max_out_files, sam_prefix, bwa_P);
#ifdef USE_MPI
	char hostname[20];
	MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	gethostname(hostname, 20);
	printf ("Rank %d/%d on %s\n", rank, numtasks, hostname);
	if (numtasks != bwa_P)
	{
		printf ("!!! Num.Ranks should be same as #procs used for Alignment\n");
		MPI_Finalize();
		return 0;
	}
#endif

	if (rank == 0)
		printf ("NumTasks %d rank %d\n", numtasks, rank);
	char line[200];
	NUM_CONTIGS=0;
        fgets(line, 200, fp_fai);
        while (!feof(fp_fai))
        {
		NUM_CONTIGS++;
		fgets(line, 200, fp_fai);
        }
	fclose(fp_fai);

	fp_fai = fopen(fai_file, "r");
        long length=0;
	fileList = (fileInfo *)malloc(NUM_CONTIGS * sizeof(fileInfo));
        fgets(line, 200, fp_fai);
        int contigNo=0;
	char contigName[100];
	char fileName[500];
	long totalRecords=0;
	long genomeSize=0;
        while (!feof(fp_fai))
        {
		sscanf(line, "%s\t%ld", contigName, &length);
		int fLen = strlen(contigName);
		fileList[contigNo].fileName = (char *)malloc(fLen+1);
		strcpy(fileList[contigNo].fileName, contigName);
		fileList[contigNo].numRecords = 0;
		fileList[contigNo].seqLength = length;
		genomeSize += fileList[contigNo].seqLength;
		if (maxSeqLen < length)
			maxSeqLen = length;
		contigNo++;
		fgets(line, 200, fp_fai);
        }
	fclose(fp_fai);
	NUM_CONTIGS = contigNo;

	split_fasta (max_out_files, fileList, NUM_CONTIGS, &numOutChunks, &numMaxChunks, &maxChunkSize);
	if (rank == 0)
		printf ("NUM_CONTIGS %d numOutChunks %d numMaxChunks %d maxChunkSize %d\n", NUM_CONTIGS, numOutChunks, numMaxChunks, maxChunkSize);

	struct timeval stime, etime, time1, time2, time3, time4;
	gettimeofday (&stime, NULL);
	size_t *segSizes_aeb = (size_t *)malloc(NUM_CONTIGS*numMaxChunks*numtasks*sizeof(size_t));
	assert (segSizes_aeb != NULL);
	size_t *segSizes_aib = (size_t *)malloc(NUM_CONTIGS*numMaxChunks*numtasks*sizeof(size_t));
	assert (segSizes_aib != NULL);
	int segNo = 0;
	i=0;
	for (segNo=0; segNo < NUM_CONTIGS*numMaxChunks; segNo++)
	{
		sprintf (fileName, "%s%d/out_s%d.aeb", sam_prefix, rank, segNo);
		int status=stat (fileName, &buf);
		if (status == -1)
			segSizes_aeb[rank*NUM_CONTIGS*numMaxChunks+segNo] = 0;
		else
			segSizes_aeb[rank*NUM_CONTIGS*numMaxChunks+segNo] = (buf.st_size/sizeof(fullRec));

		sprintf (fileName, "%s%d/out_s%d.aib", sam_prefix, rank, segNo);
		status=stat (fileName, &buf);
		if (status == -1)
			segSizes_aib[rank*NUM_CONTIGS*numMaxChunks+segNo] = 0;
		else
			segSizes_aib[rank*NUM_CONTIGS*numMaxChunks+segNo] = (buf.st_size/sizeof(otherRec));
	}
#ifdef USE_MPI
	MPI_Allgather(MPI_IN_PLACE, NUM_CONTIGS*numMaxChunks, MPI_UNSIGNED_LONG, segSizes_aeb, NUM_CONTIGS*numMaxChunks, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
	MPI_Allgather(MPI_IN_PLACE, NUM_CONTIGS*numMaxChunks, MPI_UNSIGNED_LONG, segSizes_aib, NUM_CONTIGS*numMaxChunks, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
	int numValidSegs = 0;
	int *segIds = (int *)malloc(NUM_CONTIGS*numMaxChunks*sizeof(int));
	assert (segIds != NULL);
	size_t *totSizes_aeb = (size_t *)malloc(NUM_CONTIGS*numMaxChunks*sizeof(size_t));
	assert (totSizes_aeb != NULL);
	size_t *totSizes_aib = (size_t *)malloc(NUM_CONTIGS*numMaxChunks*sizeof(size_t));
	assert (totSizes_aib != NULL);
	for (segNo=0; segNo < NUM_CONTIGS*numMaxChunks; segNo++)
	{
		size_t aebSize=0, aibSize=0;
		for (i=0; i<numtasks; i++)
		{
			aebSize += segSizes_aeb[i*NUM_CONTIGS*numMaxChunks + segNo];
			aibSize += segSizes_aib[i*NUM_CONTIGS*numMaxChunks + segNo];
		}
		if (aebSize>0 || aibSize>0)
		{
			segIds[numValidSegs] = segNo;
			totSizes_aeb[numValidSegs] = aebSize;
			totSizes_aib[numValidSegs] = aibSize;
			numValidSegs++;
		}
	}
	for (segNo=0; segNo < numValidSegs; segNo+=numtasks)
	{
		fullRec *Fr = NULL;
		otherRec *Or = NULL;
		fullRec *Fr2 = NULL;
		otherRec *Or2 = NULL;
		long comm_time = 0;
		size_t aebWriteVol=0, aibWriteVol=0;
		gettimeofday (&time1, NULL);
		for (i=0; i<numtasks; i++) {
			if (segNo+i >= numValidSegs) break;
			if (i==rank)
			{
				Fr = (fullRec *)malloc(totSizes_aeb[segNo+i]*sizeof(fullRec));
				assert (Fr != NULL);
				Or = (otherRec *)malloc(totSizes_aib[segNo+i]*sizeof(otherRec));
				assert (Or != NULL);
				Fr2 = (fullRec *)malloc(totSizes_aeb[segNo+i]*sizeof(fullRec));
				assert (Fr2 != NULL);
				Or2 = (otherRec *)malloc(totSizes_aib[segNo+i]*sizeof(otherRec));
				assert (Or2 != NULL);
			}
			//if (segIds[segNo+i] != 105)
			//	continue;
			int *recvCounts = (int *)malloc(numtasks*sizeof(int));
			assert (recvCounts != NULL);
			int *displ = (int *)malloc(numtasks*sizeof(int));
			assert (displ != NULL);
			size_t *recvCounts2 = (size_t *)malloc(numtasks*sizeof(size_t));
			assert (recvCounts2 != NULL);
			size_t *displ2 = (size_t *)malloc(numtasks*sizeof(size_t));
			assert (displ2 != NULL);
			int j=0;
			for (j=0; j<numtasks; j++) {
				recvCounts[j] = segSizes_aeb[j*NUM_CONTIGS*numMaxChunks + segIds[segNo+i]]*sizeof(fullRec);
				recvCounts2[j] = segSizes_aeb[j*NUM_CONTIGS*numMaxChunks + segIds[segNo+i]]*sizeof(fullRec);
				if (j==0) {
					displ[j] = 0;
					displ2[j] = 0;
				}
				else {
					displ[j] = displ[j-1] + recvCounts[j-1];
					displ2[j] = displ2[j-1] + recvCounts2[j-1];
				}
			}
			fullRec *lFr = NULL;
			sprintf (fileName, "%s%d/out_s%d.aeb", sam_prefix, rank, segIds[segNo+i]);
			FILE *fp_sam = fopen(fileName, "rb");
	                if(fp_sam != NULL)
			{
				lFr = (fullRec *)malloc(recvCounts2[rank]);
				assert (lFr != NULL);
				size_t nread = fread (lFr, sizeof(fullRec), recvCounts2[rank]/sizeof(fullRec), fp_sam);
				//printf ("rank %d fname %s exp_size %lu displ %lu act %lu\n", rank, fileName, recvCounts2[rank]/sizeof(fullRec), displ2[rank]/sizeof(fullRec), nread);
				assert (nread == recvCounts2[rank]/sizeof(fullRec));
				fclose (fp_sam);
			}
			//MPI_Gatherv(lFr, recvCounts[rank], MPI_UNSIGNED_CHAR, Fr, recvCounts, displ, MPI_UNSIGNED_CHAR, i, MPI_COMM_WORLD);
			mpi_gatherv(lFr, recvCounts2[rank], Fr2, recvCounts2, displ2, i, numtasks, rank);

			for (j=0; j<numtasks; j++) {
				recvCounts[j] = segSizes_aib[j*NUM_CONTIGS*numMaxChunks + segIds[segNo+i]]*sizeof(otherRec);
				recvCounts2[j] = segSizes_aib[j*NUM_CONTIGS*numMaxChunks + segIds[segNo+i]]*sizeof(otherRec);
				if (j==0) {
					displ[j] = 0;
					displ2[j] = 0;
				}
				else {
					displ[j] = displ[j-1] + recvCounts[j-1];
					displ2[j] = displ2[j-1] + recvCounts2[j-1];
				}
			}
			sprintf (fileName, "%s%d/out_s%d.aib", sam_prefix, rank, segIds[segNo+i]);
			otherRec *lOr = NULL;
			fp_sam = fopen(fileName, "rb");
	                if(fp_sam != NULL)
			{
				lOr = (otherRec *)malloc(recvCounts2[rank]);
				assert (lOr != NULL);
				assert (fread (lOr, sizeof(otherRec), recvCounts2[rank]/sizeof(otherRec), fp_sam) == recvCounts2[rank]/sizeof(otherRec));
				fclose (fp_sam);
				//int ret=unlink(fileName);
	                        //if (ret==-1)
	                        //        printf ("File %s Not Deleted\n", fileName);
			}
			//MPI_Gatherv(lOr, recvCounts[rank], MPI_UNSIGNED_CHAR, Or, recvCounts, displ, MPI_UNSIGNED_CHAR, i, MPI_COMM_WORLD);
			comm_time += mpi_gatherv(lOr, recvCounts2[rank], Or2, recvCounts2, displ2, i, numtasks, rank);
			// Process seg segNo+i
			free(recvCounts);
			free(displ);
			free(recvCounts2);
			free(displ2);
			if (lFr) free (lFr); lFr = NULL;
			if (lOr) free (lOr); lOr = NULL;

		}
		gettimeofday (&time2, NULL);
		if (Fr2)
			qsort(Fr2, totSizes_aeb[segNo+rank], sizeof(fullRec), sizecmpfn_aeb);
		if (Or2)
			qsort(Or2, totSizes_aib[segNo+rank], sizeof(otherRec), sizecmpfn_aib);
		gettimeofday (&time3, NULL);
		if (Fr2)
		{
			int segId = segIds[segNo+rank];
			int contigNo = segId/numMaxChunks;
			segId = segId%numMaxChunks;
			sprintf (fileName, "%s0/C%d_%d_sorted.aeb", sam_prefix, contigNo, segId);
			
			FILE *fp_out = fopen(fileName, "wb");
			assert(fp_out!=NULL);
			fwrite (Fr2, sizeof(fullRec), totSizes_aeb[segNo+rank], fp_out);
			fclose (fp_out);
			aebWriteVol = totSizes_aeb[segNo+rank];
		}
		if (Or2)
		{
			int segId = segIds[segNo+rank];
			int contigNo = segId/numMaxChunks;
			segId = segId%numMaxChunks;
			sprintf (fileName, "%s0/C%d_%d_sorted.aib", sam_prefix, contigNo, segId);
	
			FILE *fp_out = fopen(fileName, "wb");
			assert(fp_out!=NULL);
			fwrite (Or2, sizeof(otherRec), totSizes_aib[segNo+rank], fp_out);
			fclose (fp_out);
			aibWriteVol = totSizes_aib[segNo+rank];
		}
		gettimeofday (&time4, NULL);
		if (Fr2 || Or2)
		{
			int segId = segIds[segNo+rank];
			int contigNo = segId/numMaxChunks;
			segId = segId%numMaxChunks;
			
			write_histfile(sam_prefix, segNo+rank, contigNo, segId, totSizes_aeb, totSizes_aib, Fr2, Or2);

			long readT = (time2.tv_sec * 1000000 + time2.tv_usec) - (time1.tv_sec * 1000000 + time1.tv_usec);
			long sortT =  (time3.tv_sec * 1000000 + time3.tv_usec) - (time2.tv_sec * 1000000 + time2.tv_usec);
			long writeT = (time4.tv_sec * 1000000 + time4.tv_usec) - (time3.tv_sec * 1000000 + time3.tv_usec);
			printf ("[%d] C%d_%d Read+comm(%lf sec) %lf sec , Sort %lf sec, Write %lf sec (AEB %lu AIB %lu)\n", rank, contigNo, segId, (double)comm_time/1000000, (double)readT/1000000, (double)sortT/1000000, (double)writeT/1000000, aebWriteVol*sizeof(fullRec), aibWriteVol*sizeof(otherRec));
		}
		if (Fr) free (Fr);
		if (Or) free (Or);
		Fr = NULL; Or = NULL;
		if (Fr2) free (Fr2);
		if (Or2) free (Or2);
		Fr2 = NULL; Or2 = NULL;
	}
	free (segSizes_aeb);
	free (segSizes_aib);
	free(segIds);
	free (totSizes_aeb);
	free (totSizes_aib);

	gettimeofday (&etime, NULL);
	long elapsed = (etime.tv_sec * 1000000 + etime.tv_usec) - (stime.tv_sec * 1000000 + stime.tv_usec);
	printf ("[%d] - %lf sec\n", rank, (double)elapsed/1000000);
//	free (nrecs);
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	for (i=0; i<NUM_CONTIGS; i++)
		free(fileList[i].fileName);
	free(fileList);
#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
