#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int dim=2;

/*** fills the string with a random genomic sequence ***/

void populate (char* string, long long int length)
{
	int random;
	long long int i;
	for(i=0;i<length-1;i++)
	{
		random=rand()%4;
		switch(random)
		{
			case 0: string[i]='A'; break;
			case 1: string[i]='C'; break;
			case 2: string[i]='T'; break;
			case 3: string[i]='G';
		}
		if (i%100000==0)
			srand(rand()); // re-randomize the seed
	}
	string[length-1]='\0';
}

/*** function which concretely finds all the correspondences, returning their index ***/

long long int* matcher (char* sequence, char* pattern, int* count, int hole, int slen, int plen)
{
	long long int i, j, k; //counters for the scan of sequence and of pattern, forward and backward
	int h=0, htot=0; //counters for the length of the current hole and of the total n of holes encountered in the current comparison
	long long int *corresp=malloc(dim*sizeof(long long int)); // lists the indices of correspondances, initially from 0 to dim
	for (i=0; i<=slen-plen; i++) // scan the sequence
	{
		h=0;
		htot=0;
		for (j=0; j<plen-1; j++) // scan the pattern forward. plen-1 because the last char is '\0'
		{
			if (sequence[i+j]!=pattern[j])
			{
				h++;
				htot++;
				if (h>hole/100 || htot>hole) // the largest contiguous hole is hole/100
					break; // abandon this comparison
			}
		}
		h=0;
		htot=0;
		for (k=0; k<plen-1; k++) // scan the pattern backward. plen-1 because the last char is '\0'
		{
			
			if (sequence[i+k]!=pattern[plen-2-k])
			{
				h++;
				htot++;
				if (h>hole/100 || htot>hole)
					break; // abandon this comparison
			}
		}
		if (j==plen-1||k==plen-1) // all the pattern matches (at least one of the two)
		{
			if (*count>=dim) // max of corrispondences reached...
			{
				dim*=2;
				corresp=(long long int*)realloc (corresp, dim*sizeof(long long int)); // ...reallocate with doubled dim
				if (!corresp)
				{
					perror("realloc failed");
					//return NULL;
				}
			}
			corresp[*count]=i; // add the new correspondence index
			(*count)++; // number of correspondences found
		}
	}
	corresp=(long long int*)realloc(corresp, *count*sizeof(long long int)); // trim allocated memory
	return corresp;
}

/*** main ***/

int main(int argc, char *argv[])
{
	double totaltime=omp_get_wtime();
	srand(time(NULL));
	int i=0, count=0;
	long long int slen=10000001, plen=10001; // 64bit variables
	double seconds;
	if (argc==3)
	{
		slen=atoi(argv[1]); // length of sequence
		plen=atoi(argv[2]); // length of pattern
		if (plen>slen)
		{
			perror("inconsistent data");
			return -1;
		}
	}
	long long int *corresp;
	char *sequence=malloc(slen*sizeof(char)); // sequence to search in
	char *pattern=malloc(plen*sizeof(char)); // substring to be found
	if(!sequence || !pattern)
	{
		perror("malloc failed");
		return -1;
	}
	seconds=omp_get_wtime();
	populate (sequence, slen);
	populate (pattern, plen);
	//printf("sequence: %s\n",sequence);
	//printf("pattern: %s\n",pattern);
	seconds=omp_get_wtime()-seconds;
	printf("Peopling time: %lf\n", seconds);
	
	seconds=omp_get_wtime();
	corresp=matcher(sequence, pattern, &count, plen/10, slen, plen);
	seconds=omp_get_wtime()-seconds;
	
	printf("%d forward and backward correspondences found\n", count);
	puts("Positions: ");
	for(i=0;i<count;i++)
		printf("%lli  ", corresp[i]);
	printf("\nComputation time: %lf\n", seconds);
	free(sequence);
	free(pattern);
	free(corresp);
	printf("Total time: %lf\nNumber of cores: %d\n", omp_get_wtime()-totaltime, omp_get_num_threads());
	return 0;
}