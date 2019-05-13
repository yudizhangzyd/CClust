#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "align.h"
#include "error.h"

/**
 * Needleman-Wunsch alignment.
 *
 * Note: input sequences must end with string termination character, '\0'
 *
 * @param s1	first sequence
 * @param s2	second sequence
 * @param len1	length of first sequence
 * @param len2	length of first sequence
 * @param score	scores
 * @param gap_p
 * @param band	band, -1 for no band
 * @param ends_free	ends-free alignment
 * @param perr	error probability in read (second sequence)
 * @param alen	pointer to alignment length
 * @return	alignment
 */
unsigned char **nwalign(unsigned char const * const s1, unsigned char const * const s2,
	size_t len1, size_t len2, int score[4][4], int gap_p, int band,
	int ends_free, double const *perr, int *err, size_t *alen) {
	static size_t nnw = 0;
	size_t i, j;
	size_t l, r;
	size_t iband = band >= 0 ? band : 0;
	double diag, left, up;
/*
	unsigned int len1 = strlen(s1);
	unsigned int len2 = strlen(s2);
*/

	*err = NO_ERROR;

	unsigned int nrow = len1 + 1;
	unsigned int ncol = len2 + 1;

	//int *d = (int *) malloc(nrow * ncol * sizeof(int)); //E
	double *d = (double *) malloc(nrow * ncol * sizeof(double)); //E
	int *p = (int *) malloc(nrow * ncol * sizeof(int)); //E
	if (d == NULL || p == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "d & p");
		return NULL;
	}

	// Fill out left columns of d, p.
	for (i = 0; i <= len1; i++) {
		d[i*ncol] = ends_free ? 0 : i * gap_p; // ends-free gap
		p[i*ncol] = 3;
	}

	// Fill out top rows of d, p.
	for (j = 0; j <= len2; j++) {
		d[j] = ends_free ? 0 : j * gap_p; // ends-free gap
		p[j] = 2;
	}

	// Calculate left/right-bands in case of different lengths
	size_t lband, rband;
	if (len2 > len1) {
		lband = iband;
		rband = iband + len2 - len1;
	} else if (len1 > len2) {
		lband = iband + len1 - len2;
		rband = iband;
	} else {
		lband = iband;
		rband = iband;
	}

	// Fill out band boundaries of d.
	if (band >= 0 && (iband < len1 || iband < len2)) {
		for (i = 0; i <= len1; i++) {
			if ((int) i - (int) lband - 1 >= 0)
				d[i*ncol + i - lband - 1] = -9999;
			if (i + rband + 1 <= len2)
				d[i*ncol + i + rband + 1] = -9999;
		}
	}

	// Fill out the body of the DP matrix.
	for (i = 1; i <= len1; i++) {
		if (band >= 0) {
			l = i - lband;
			if (l < 1)
				l = 1;
			r = i + rband;
			if (r > len2)
				r = len2;
		} else {
			l = 1;
			r = len2;
		}

		for (j = l; j <= r; j++) {
			// Score for the left move.
			if (i == len1)
				left = d[i*ncol + j - 1]
					+ (ends_free ? 0 : gap_p); // Ends-free gap.
			else
				left = d[i*ncol + j - 1] + gap_p;

			// Score for the up move.
			if (j == len2)
				up = d[(i-1)*ncol + j]
					+ (ends_free ? 0 : gap_p); // Ends-free gap.
			else
				up = d[(i-1)*ncol + j] + gap_p;

			// Score for the diagonal move.
			diag = d[(i-1)*ncol + j-1]
				+ (perr ? perr[j-1] : 1.) * score[(int) s1[i-1]][(int) s2[j-1]];

			// Break ties and fill in d,p.
			if (up >= diag && up >= left) {
				d[i*ncol + j] = up;
				p[i*ncol + j] = 3;
			} else if (left >= diag) {
				d[i*ncol + j] = left;
				p[i*ncol + j] = 2;
			} else {
				d[i*ncol + j] = diag;
				p[i*ncol + j] = 1;
			}
		}
	}

	unsigned char *al0 = (unsigned char *) malloc((len1+len2) * sizeof(unsigned char));
	unsigned char *al1 = (unsigned char *) malloc((len1+len2) * sizeof(unsigned char));
	if (al0 == NULL || al1 == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al0 & al1");
		return NULL;
	}

	// Trace back over p to form the alignment.
	size_t len_al = 0;
	i = len1;
	j = len2;

	while ( i > 0 || j > 0 ) {
		switch ( p[i*ncol + j] ) {
			case 1:
				al0[len_al] = s1[--i];
				al1[len_al] = s2[--j];
				break;
			case 2:
				al0[len_al] = '-';
				al1[len_al] = s2[--j];
				break;
			case 3:
				al0[len_al] = s1[--i];
				al1[len_al] = '-';
				break;
			default:
				*err = OUT_OF_BAND;
				mmessage(ERROR_MSG, *err,
					"NW alignment out of range");
				return NULL;
		}
		len_al++;
	}
//	al0[len_al] = '\0';
//	al1[len_al] = '\0';
/*
	for (size_t i = 0; i < nrow; ++i) {
		for (size_t j = 0; j < ncol; ++j)
			fprintf(stderr, " %2d", d[i*ncol + j]);
		fprintf(stderr, "\n");
	}
*/


	// Allocate memory to alignment strings.
	unsigned char **al = (unsigned char **) malloc( 2 * sizeof(unsigned char *) ); //E
	if (al == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al");
		return NULL;
	}

	al[0] = (unsigned char *) malloc(len_al); //E
	al[1] = (unsigned char *) malloc(len_al); //E
	if (al[0] == NULL || al[1] == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al[]");
		return NULL;
	}

	// Reverse the alignment strings (since traced backwards).
	for (i = 0 ; i < len_al ; i++) {
		al[0][i] = al0[len_al-i-1];
		al[1][i] = al1[len_al-i-1];
	}
//	al[0][len_al] = '\0';
//	al[1][len_al] = '\0';

	// Free allocated memory
	free(d);
	free(p);
	free(al0);
	free(al1);

	*alen = len_al;

	nnw++;
	return al;
} /* nwalign */
