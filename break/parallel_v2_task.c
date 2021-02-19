#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define H_ERR 10
#define H_WIDTH 1000

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
	long long int i = 0, j = 0, k = 0; // counters: 'j' for fwd, 'k' for backwd
	// counters for width of current hole and number of holes
	int count_L_hole = 0, err_cnt = 0;
	int dim = 2;

	long long int *corresp = malloc(dim * sizeof(long long int));
	long long int max_holes = plen / H_ERR; // H_ERR = 10
	long long int max_L_hole = plen / H_WIDTH; // H_WIDTH = 1000

	#pragma omp parallel
	{
	#pragma omp single
	{
	#pragma omp task firstprivate(count_L_hole,err_cnt) private(i,j)
	{
	#pragma omp parallel for private (count_L_hole,err_cnt,j) schedule (dynamic,plen)

	for (i=0; i<=slen-plen; i++) {

		if (i == 0) {
			printf("NUM_THREADS %i\n", omp_get_num_threads());
			printf("Chunk dimension %i\n", plen / omp_get_num_threads());
		}

		count_L_hole = 0;
		err_cnt = 0;

		// Forward:
		for (j = 0; j < plen-1; j++) {
			if (sequence[i+j] != pattern[j]) {
				count_L_hole++;

				if (count_L_hole == 1) // starting the hole...
					err_cnt++;

				if (count_L_hole > max_L_hole || err_cnt > max_holes)
					break;
			} else {
				count_L_hole = 0;
			}
		}

		count_L_hole = 0;
		err_cnt = 0;

		// Backward:
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
			#pragma omp critical(realloc)
			{
			if (*count >= dim) {
				dim *= 2;
				corresp=(long long int*) realloc(corresp, dim*sizeof(long long int));
				if (!corresp) {
					fprintf(stderr, "Realloc failed");
				}
			}
			} // end Critical realloc

			#pragma omp critical(update)
			{
			corresp[*count] = i;
			(*count)++;
			} // end Critical update
		}
	} /* end for loop */
	} /* end pragma task */
	} /* end pragma single */
	} /* end pragma parallel */

	corresp = (long long int*) realloc(corresp, *count * sizeof(long long int)); // trim allocated memory
	return corresp;
}

int main() {
	long long int seq_len, pat_len;
	long long int *corresp;
 	int count = 0;

	int n_thr;

	// User-acquired inputs:
	printf("\n* Insert desired sequence length: ");
	scanf("%lli", &seq_len);
	printf("\n* Insert desired pattern length: ");
	scanf("%lli", &pat_len);
	printf("\n* Insert number of threads: ");
	scanf("%i", &n_thr);

	omp_set_num_threads(n_thr);

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

	double init_populate = omp_get_wtime();

	// (Random) generation of sequence and pattern
	populate(sequence, seq_len);
	populate(pattern, pat_len);

	double end_populate = omp_get_wtime() - init_populate;
	double init_match = omp_get_wtime();

 	// Fwd sequencing with holes:
	corresp = hole_forward(sequence, pattern, &count, seq_len, pat_len);

	double end_match = omp_get_wtime() - init_match;

	printf("\n%i correspondances found!\n", count);

	printf("\nPositions of the matches:\n");
	for (int i = 0; i < count; i++) {
		printf("* %lli\n", corresp[i] + 1);
	}

	double total_time_end = omp_get_wtime() - total_time_init;

  	printf("\n* POPULATION TIME: %.5f\n\n", end_populate);
  	printf("\n* MATCH TIME: %.5f\n\n", end_match);
  	printf("\n* TOTAL TIME: %.5f\n\n", total_time_end);

	free(sequence);
	free(pattern);
	free(corresp);

	return 0;
}
