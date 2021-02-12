#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define H_ERR 10
#define H_WIDTH 1000

// Struct 'conteggio' contains number of matches as
// well as their position in the sequence
typedef struct {
	int count_fwd, count_bck;
	int tot_count;
	long long int *corresp;
} conteggio;

// Randomly generate a sequence of genomes:
void populate(char str[], long long int length) {
	const char az_base[4] = {'A', 'C', 'G', 'T'};
	// Remove after testing
	srand(4);

	for (int i=0; i<length; i++) {
		int r = rand() % 4;
		str[i] = az_base[r];
	}

	str[length] = '\0';
}

// Fwd & Bckwd sequencing in case of genes' undesired modifications
void sequencing(char seq[], char pat[], long long int s_len, long long int p_len, conteggio *c) {
	// 'i' index: sequence starting position;
	// 'j' index: pattern vs sub-sequence comparison.
	int i = 0, j = 0, dim = 2;
	int err_cnt = 0, count_L_hole = 0;
	int max_holes = (int) p_len / H_ERR;
	int max_L_hole = (int) p_len / H_WIDTH;

	c->count_fwd = c->count_bck = c->tot_count = 0;
	c->corresp = malloc(dim * sizeof(long long int));

	if (!c->corresp) {
		fprintf(stderr, "Malloc failed\n");
	}

	// Forward: //
	for (i = 0; i <= s_len-p_len; i++) {
		for (j = 0; j < p_len; j++) {

			// Realloc first, if needed
			if ((c->tot_count) >= dim) {
				dim *= 2;
				c->corresp = realloc(c->corresp, dim * sizeof(long long int));

				if (!c->corresp) {
					fprintf(stderr, "Realloc failed!\n");
				}
			}

			// Count only if the comparison is true given the
			// max number of allowed 'holes':
			if (err_cnt <= max_holes && count_L_hole <= max_L_hole) {
				if (seq[i+j] != pat[j]) {
					count_L_hole++;

					if (count_L_hole == 1)
						err_cnt++;

					if (count_L_hole > max_L_hole) {
						break;
					}
				} else {
					count_L_hole = 0;
				}
				
				if (j == p_len-1 && err_cnt <= max_holes) {
					c->corresp[c->tot_count] = i;
					printf("\n** DEBUGGONE:\n");
					printf("\nCounter - Position: %i - %lli\n", c->tot_count, c->corresp[c->tot_count]);
					(c->count_fwd)++;
					(c->tot_count)++;
				}
			}
		}

		err_cnt = 0;
		count_L_hole = 0;
	}

	// Backward: //
	for (i = s_len-1; i >= p_len-1; i--) {
		for (j = p_len-1; j >= 0; j--) {

			// Realloc first, if needed
			if ((c->tot_count) >= dim) {
				dim *= 2;
				c->corresp = realloc(c->corresp, dim * sizeof(long long int));

				if (!c->corresp) {
					fprintf(stderr, "Realloc failed!\n");
				}
			}

			// Count only if the comparison is true given the
			// max number of allowed 'holes':
			if (err_cnt <= max_holes && count_L_hole <= max_L_hole) {
				if (seq[i-j] != pat[j]) {
					count_L_hole++;

					if (count_L_hole == 1)
						err_cnt++;

					if (count_L_hole > max_L_hole) {
						break;
					}
				} else {
					count_L_hole = 0;
				}
				
				// No violation on holes constraints: correspondance found!
				if (j == 0 && err_cnt <= max_holes) {
					c->corresp[c->tot_count] = i;
					printf("\n** DEBUGGONE:\n");
					printf("\nCounter - Position: %i - %lli\n", c->tot_count, c->corresp[c->tot_count]);
					(c->count_bck)++;
					(c->tot_count)++;
				}
			}
		}

		err_cnt = 0;
		count_L_hole = 0;
	}

	// Trim memory:
	c->corresp = realloc(c->corresp, c->tot_count * sizeof(long long int));
}

int main() {

	long long int seq_len, pat_len;
	conteggio c;
	
	// User-acquired inputs:
	printf("\n* Insert desired sequence length: ");
	scanf("%lli", &seq_len);
	printf("\n* Insert desired pattern length: ");
	scanf("%lli", &pat_len);

	/* Start timer */
	double total_time_init = omp_get_wtime();

	printf("\n==> Maximum allowed number of holes: %i\n", (int) pat_len / H_ERR);
	
	char *sequence = malloc(seq_len * sizeof(char));
	char *pattern = malloc(pat_len * sizeof(char));

	if(!sequence || !pattern) {
		fprintf(stderr, "Malloc failed\n");
		return 2;
	}

	if (pat_len > seq_len) {
		fprintf(stderr, "ERROR ON SELECTED SIZES!\n");
		return 4;
	}

	double init_populate = omp_get_wtime();

	// (Random) generation of sequence and pattern
	populate(sequence, seq_len);
	populate(pattern, pat_len);

	double end_populate = omp_get_wtime();

	// uncomment to check patterns for errors
	//printf("\n\nSEQUENCE: %s;\nPATTERN: %s\n\n", sequence, pattern);
 
	double init_match = omp_get_wtime();
	
 	// Fwd & Bckwd sequencing with holes:
	sequencing(sequence, pattern, seq_len, pat_len, &c);
	double end_match = omp_get_wtime();

	printf("\n%i correspondances found! (%i forward, %i backward)\n", c.tot_count, c.count_fwd, c.count_bck);

	printf("\nPositions of the matches:\n");
	for (int i = 0; i < c.tot_count; i++) {
		printf("* %lli\n", c.corresp[i] + 1);
	}

	double total_time_end = omp_get_wtime();

  	printf("\nPOPULATION TIME: %.5f\n\n", end_populate - init_populate);
  	printf("\nMATCH TIME: %.5f\n\n", end_match - init_match);
  	printf("\nTOTAL TIME: %.5f\n\n", total_time_end - total_time_init);

	free(sequence);
	free(pattern);
	free(c.corresp);

	return 0;
}