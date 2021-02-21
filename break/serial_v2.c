#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define H_ERR 10
#define H_WIDTH 1000

int dim = 2;

// Randomly generate a sequence of genomes:
void generate(char str[], long long int length) {
	const char az_base[4] = {'A', 'C', 'G', 'T'};
	// srand(4);

	for (int i=0; i<length; i++) {
		int r = rand() % 4;
		str[i] = az_base[r];
	}

	str[length] = '\0';
}

long long int* sequencing(char* sequence, char* pattern, int* count, int slen, int plen) {
	long long int i = 0, j = 0, k = 0; // counters: 'j' for fwd, 'k' for backwd
	int count_L_hole = 0, err_cnt = 0;

	long long int *corresp = malloc(dim * sizeof(long long int));
	long long int max_holes = plen / H_ERR; // H_ERR = 10
	long long int max_L_hole = plen / H_WIDTH; // H_WIDTH = 1000
	
	for (i = 0; i <= slen-plen; i++) {

		count_L_hole = 0;
		err_cnt = 0;

		/* Forward */
		for (j = 0; j < plen-1; j++) {
			if (sequence[i+j] != pattern[j]) {
				count_L_hole++;

				if (count_L_hole == 1) // starting the hole...
					err_cnt++;

				if (count_L_hole > max_L_hole || err_cnt > max_holes)
					break; // abandon comparison
			} else {
				count_L_hole = 0;
			}
		}

		count_L_hole = 0;
		err_cnt = 0;

		/* Backward */
		for (k = 0; k < plen-1; k++) {
			if (sequence[i+k] != pattern[plen-2-k]) {
				count_L_hole++;

				if (count_L_hole == 1)
					err_cnt++;

				if (count_L_hole > max_L_hole || err_cnt > max_holes)
					break;
			} else {
				count_L_hole = 0;
			}
		}
		
		if (j == plen-1 || k == plen-1) // new corresp. found if we get to end of pattern!
		{
			if (*count >= dim) {
				dim *= 2;
				corresp=(long long int*)realloc (corresp, dim*sizeof(long long int));
				if (!corresp) {
					perror("realloc failed");
				}
			}

			corresp[*count]=i;
			(*count)++;
		}
	}

	corresp = (long long int*) realloc(corresp, *count * sizeof(long long int)); // trim allocated memory
	return corresp;
}

void print_speedup(double f, int n) {
	double s = 0;

	s = n + (1 - n) * f;
	s = n / s;
  	printf("Expected speedup for %i cores: %.5f\n\n", n, s);
}

int main() {
	long long int seq_len, pat_len;
	long long int *corresp;
 	int count = 0;
	double f = 0;

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

	if(!sequence || !pattern) {
		fprintf(stderr, "Malloc failed\n");
		return -2;
	}

	if (pat_len > seq_len) {
		fprintf(stderr, "ERROR ON SELECTED SIZES!\n");
		return -4;
	}

	double init_generate = omp_get_wtime();

	// (Random) generation of sequence and pattern
	generate(sequence, seq_len);
	generate(pattern, pat_len);

	double end_generate = omp_get_wtime() - init_generate;

	double init_match = omp_get_wtime();

 	// Fwd sequencing with holes:
	corresp = sequencing(sequence, pattern, &count, seq_len, pat_len);

	double end_match = omp_get_wtime() - init_match;

	printf("\n%i correspondances found!\n", count);

	printf("\nPositions of the matches:\n");
	for (int i = 0; i < count; i++) {
		printf("* %lli\n", corresp[i] + 1);
	}

	double total_time_end = omp_get_wtime() - total_time_init;

  	printf("\nPOPULATION TIME: %.5f\n\n", end_generate);
  	printf("\nMATCH TIME: %.5f\n\n", end_match);
  	printf("\nTOTAL TIME: %.5f\n\n", total_time_end);

  	f = end_match / total_time_end; // Amdahl law fraction

  	printf("* Percentage: %.5f\n\n", f * 100);

  	// Expected speedups
  	int n[5] = {2, 4, 8, 16, 24};

	for (int i = 0; i < 5; i++) {
		print_speedup(f, n[i]);
	}

	free(sequence);
	free(pattern);
	free(corresp);
	
	return 0;
}
