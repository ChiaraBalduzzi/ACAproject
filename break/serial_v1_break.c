#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define H_ERR 10
#define H_WIDTH 1000

int dim = 2;

// Struct 'conteggio' contains number of matches as
// well as their position in the sequence
typedef struct {
	int count;
	long long int *corresp;
} conteggio;

// Randomly generate a sequence of genomes:
void populate(char str[], long long int length) {
	const char az_base[4] = {'A', 'C', 'G', 'T'};
	// srand(4);

	for (int i=0; i<length; i++) {
		int r = rand() % 4;
		str[i] = az_base[r];
	}

	str[length] = '\0';
}

long long int* hole_forward(char* sequence, char* pattern, int* count, int slen, int plen) {
	long long int i = 0, j = 0;
	int count_L_hole = 0, err_cnt = 0;
	long long int *corresp = malloc(dim * sizeof(long long int));
	long long int max_holes = plen / H_ERR; // H_ERR = 10
	long long int max_L_hole = plen / H_WIDTH; // H_WIDTH = 1000
	
	for (i=0; i<=slen-plen; i++) // scan the sequence
	{
		count_L_hole = 0;
		err_cnt = 0;
		for (j = 0; j < plen-1; j++) // scan the pattern forward. plen-1 because the last char is '\0'
		{
			if (sequence[i+j] != pattern[j]) {
				count_L_hole++;
				//err_cnt++;

				if (count_L_hole == 1) // starting the hole...
					err_cnt++;

				if (count_L_hole > max_L_hole || err_cnt > max_holes) // the largest contiguous hole is max_holes/100
					break; // abandon this comparison
			} else {
				count_L_hole = 0;
			}
		}
		
		if (j == plen-1) // all the pattern matches (at least one of the two)
		{
			if (*count >= dim) // max of corrispondences reached...
			{
				dim *= 2;
				corresp=(long long int*)realloc (corresp, dim*sizeof(long long int)); // ...reallocate with doubled dim
				if (!corresp) {
					perror("realloc failed");
					//return NULL;
				}
			}
			corresp[*count]=i; // add the new correspondence index
			(*count)++; // number of correspondences found
		}
	}
	corresp = (long long int*) realloc(corresp, *count * sizeof(long long int)); // trim allocated memory
	return corresp;
}

int main() {
	long long int seq_len, pat_len;
	//conteggio c_fw, c_bck;

	// User-acquired inputs:
	printf("\n* Insert desired sequence length: ");
	scanf("%lli", &seq_len);
	printf("\n* Insert desired pattern length: ");
	scanf("%lli", &pat_len);

	/* Start timer */
	double total_time_init = omp_get_wtime();

	printf("\n==> Maximum allowed number of holes: %i\n", (int) pat_len / H_ERR);

	char *sequence = malloc((seq_len+1) * sizeof(char));
	char *pattern = malloc((pat_len+1) * sizeof(char));
	char *rev_pat = malloc((pat_len+1) * sizeof(char));

	if(!sequence || !pattern || !rev_pat) {
		fprintf(stderr, "Malloc failed\n");
		return -2;
	}

	if (pat_len > seq_len) {
		fprintf(stderr, "ERROR ON SELECTED SIZES!\n");
		return -4;
	}

	double init_populate = omp_get_wtime();

	// (Random) generation of sequence and pattern
	populate(sequence, seq_len);
	populate(pattern, pat_len);

	double end_populate = omp_get_wtime();

	for (int i = 0; i < pat_len; i++) {
		rev_pat[i] = pattern[pat_len-i-1];
	}
	rev_pat[pat_len] = '\0';

	double init_match = omp_get_wtime();

 	// Fwd sequencing with holes:
 	long long int *corresp_fw;
 	long long int *corresp_bk;
 	int count_fw = 0, count_bk = 0;
	corresp_fw = hole_forward(sequence, pattern, &count_fw, seq_len, pat_len);
	corresp_bk = hole_forward(sequence, rev_pat, &count_bk, seq_len, pat_len);

	double end_match = omp_get_wtime() - init_match;

	printf("\n%i correspondances found! (%i forward, %i backward)\n", (count_fw + count_bk), count_fw, count_bk);

	printf("\nPositions of the matches:\n");
	for (int i = 0; i < count_fw; i++) {
		printf("* %lli\n", corresp_fw[i] + 1);
	}

	// In bckwd add pat_len as corresp. idx points to the end of the comparison
	for (int i = 0; i < count_bk; i++) {
		printf("* %lli\n", corresp_bk[i] + pat_len);
	}

	double total_time_end = omp_get_wtime() - total_time_init;

  	printf("\nPOPULATION TIME: %.5f\n\n", end_populate - init_populate);
  	printf("\nMATCH TIME: %.5f\n\n", end_match);
  	printf("\nTOTAL TIME: %.5f\n\n", total_time_end);

  	printf("Percentage: %.5f\n", (end_match / total_time_end) * 100);

  	int n = 16;
  	double s = n + (1 - n) * (end_match / total_time_end);
  	s = n / s;
  	printf("Expected speedup for %i cores: %.5f\n\n", n, s);

	free(sequence);
	free(pattern);
	free(rev_pat);
	/*free(c_fw.corresp);
	free(c_bck.corresp);*/

	return 0;
}
