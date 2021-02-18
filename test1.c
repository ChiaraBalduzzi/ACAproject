#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

int dim = 2;

// Randomly generate a sequence of genomes:
void populate(char str[], long long int length) {
	const char az_base[4] = {'A', 'C', 'G', 'T'};
	// Remove after testing
	//srand(4);

	for (int i=0; i<length; i++) {
		int r = rand() % 4;
		str[i] = az_base[r];
	}

	str[length] = '\0';
}

/*** function which concretely finds all the correspondences, returning their index ***/

long long int* matcher (char* sequence, char* pattern, int* count, int hole, int slen, int plen)
{
	long long int i, j, k; //counters for the scan of sequence and of pattern, forward and backward
	int h=0, htot=0; //counters for the length of the current hole and of the total n of holes encountered in the current comparison
	long long int *corresp=malloc(dim*sizeof(long long int)); // lists the indices of correspondances, initially from 0 to dim
	
	#pragma omp parallel for schedule (dynamic, plen/omp_get_num_threads()) private (j,k,h,htot) num_threads(4)
	for (i=0; i<=slen-plen; i++) // scan the sequence
	{
		h=0;
		htot=0;
		for (j=0; j<plen; j++) // scan the pattern forward. plen-1 because the last char is '\0'
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
		for (k=0; k<plen; k++) // scan the pattern backward. plen-1 because the last char is '\0'
		{
			
			if (sequence[i+k]!=pattern[plen-1-k])
			{
				h++;
				htot++;
				if (h>hole/100 || htot>hole)
					break; // abandon this comparison
			}
		}
		if (j==plen||k==plen) // all the pattern matches (at least one of the two)
		{
			#pragma omp critical (realloc)
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
			}
			#pragma omp critical (update)
			{
			corresp[*count]=i; // add the new correspondence index
			(*count)++; // number of correspondences found
			}
		}
	}
	corresp=(long long int*)realloc(corresp, *count*sizeof(long long int)); // trim allocated memory
	printf("Number of cores: %d\n",omp_get_num_threads());
	return corresp;
}

/*** main ***/

int main(int argc, char *argv[])
{
	srand(time(NULL));
	int i = 0, count = 0;
	long long int seq_len, pat_len;
	//omp_set_num_threads(4);
	double seconds;

	// User-acquired inputs:
	printf("\n* Insert desired sequence length: ");
	scanf("%lli", &seq_len);
	printf("\n* Insert desired pattern length: ");
	scanf("%lli", &pat_len);

	double totaltime = omp_get_wtime();

	long long int *corresp;
	char *sequence=malloc((seq_len+1)*sizeof(char)); // sequence to search in
	char *pattern=malloc((pat_len+1)*sizeof(char)); // substring to be found
	if(!sequence || !pattern)
	{
		perror("malloc failed");
		return -1;
	}
	seconds=omp_get_wtime();
	populate (sequence, seq_len);
	populate (pattern, pat_len);
	//printf("sequence: %s\n",sequence);
	//printf("pattern: %s\n",pattern);
	seconds=omp_get_wtime()-seconds;
	printf("Peopling time: %lf\n", seconds);
	
	seconds=omp_get_wtime();
	corresp=matcher(sequence, pattern, &count, pat_len/10, seq_len, pat_len);
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