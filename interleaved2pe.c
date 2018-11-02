#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

int main (int argc, char **argv)
{
	if (argc != 2)
	{
		printf ("Usage: interleaved2pe <in.fq>\n");
		return 1;
	}
	int pos=strlen(argv[1])-1;
	char *suffix = &argv[1][pos];
	while ((suffix != argv[1]) && (*suffix!='.')) {suffix--; pos--;}
	if ((suffix==argv[1]) || (strcmp (suffix, ".fastq") && strcmp (suffix, ".fq")))
	{
		printf ("Invalid FASTQ file: %s\n", argv[1]);
		printf ("Expected file extension .fq or .fastq\n");
		return 1;
	}
	suffix++;
	char *outprefix = (char *)malloc (strlen(argv[1]));
	strcpy (outprefix, argv[1]);
	outprefix[pos] = '\0';
	printf ("prefix: %s\n", outprefix);

	FILE *fp_in = fopen (argv[1], "r");
	if (fp_in == NULL)
	{
		printf ("Error: FASTQ file %s does not exist\n", argv[1]);
		return 1;
	}
	char *outfile = (char *)malloc (strlen(argv[1])+2);
	sprintf (outfile, "%s.1.%s", outprefix, suffix);
	printf ("outfile: %s\n", outfile);
	FILE *fp_out[2] = {NULL, NULL};
	fp_out[0] = fopen (outfile, "w");
	assert (fp_out[0] != NULL);
	sprintf (outfile, "%s.2.%s", outprefix, suffix);
	printf ("outfile: %s\n", outfile);
	fp_out[1] = fopen (outfile, "w");
	assert (fp_out[1] != NULL);

	char line [1000];
	fgets (line, 500, fp_in);
	long readNo=0;
	while (!feof (fp_in))
	{
		fputs (line, fp_out[readNo&1]);
		fgets (line, 500, fp_in);
		assert (!feof (fp_in));
		fputs (line, fp_out[readNo&1]);
		fgets (line, 500, fp_in);
		assert (!feof (fp_in));
		fputs (line, fp_out[readNo&1]);
		fgets (line, 500, fp_in);
		assert (!feof (fp_in));
		fputs (line, fp_out[readNo&1]);
		fgets (line, 500, fp_in);
		readNo++;
	}

	fclose (fp_out[0]);
	fclose (fp_out[1]);
	fclose (fp_in);
	return 0;
}
