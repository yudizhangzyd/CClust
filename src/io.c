#include "io.h"
#include "math.h"
#include "error.h"

void fprint_doubles(FILE *fp, double *v, size_t k, int precision, int newline)
{
	size_t i;
	for (i = 0; i < k; ++i)
		fprintf(fp, " %.*f", precision, v[i]);
	if (newline) fprintf(fp, "\n");
} /* fprint_doubles */

void wprint_doubles(WINDOW *wp, double *v, size_t k, int precision, int newline)
{
	size_t i;
	for (i = 0; i < k; ++i)
		wprintw(wp, " %.*f", precision, v[i]);
	if (newline) wprintw(wp, "\n");
} /* wprint_doubles */

void fprint_ints(FILE *fp, int *v, size_t n, int width, int newline)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (width)
			fprintf(fp, " %*d", width, v[i]);
		else
			fprintf(fp, " %d", v[i]);
	if (newline) fprintf(fp, "\n");
} /* fprint_ints */

void wprint_ints(WINDOW *wp, int *v, size_t n, int width, int newline)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (width)
			wprintw(wp, " %*d", width, v[i]);
		else
			wprintw(wp, " %d", v[i]);
	if (newline) wprintw(wp, "\n");
} /* wprint_ints */

void fprint_uints(FILE *fp, unsigned int *v, size_t n, int width, int newline)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (width)
			fprintf(fp, " %*u", width, v[i]);
		else
			fprintf(fp, " %u", v[i]);
	if (newline) fprintf(fp, "\n");
} /* fprint_uints */

void wprint_uints(WINDOW *wp, unsigned int *v, size_t n, int width, int newline)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (width)
			wprintw(wp, " %*u", width, v[i]);
		else
			wprintw(wp, " %u", v[i]);
	if (newline) wprintw(wp, "\n");
} /* wprint_uints */

void fprint_size_ts(FILE *fp, size_t *v, size_t n, int width, int newline)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (width)
			fprintf(fp, " %*zu", width, v[i]);
		else
			fprintf(fp, " %zu", v[i]);
	if (newline) fprintf(fp, "\n");
} /* fprint_size_ts */

void wprint_size_ts(WINDOW *wp, size_t *v, size_t n, int width, int newline)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (width)
			wprintw(wp, " %*u", width, v[i]);
		else
			wprintw(wp, " %u", v[i]);
	if (newline) wprintw(wp, "\n");
} /* wprint_size_ts */


/**
 * Print a square nxn matrix that is stored in vectorized format to
 * file handle.
 *
 * @param fp	FILE pointer
 * @param mat	vectorized matrix
 * @param n	dimension
 * @param row	row or column vectorization
 */
void fprint_vectorized_sq_matrix(FILE *fp, double *mat, size_t n, int row)
{
	size_t i, j;

	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if ((row ? mat[i*n + j] : mat[j*n + i]) < 1e-2)
				fprintf(fp, " %8.2e",
					row ? mat[i*n + j] : mat[j*n + i]);
			else
				fprintf(fp, " %8.3f",
					row ? mat[i*n + j] : mat[j*n + i]);
		}
		fprintf(fp, "\n");
	}
} /* fprint_vectorized_sq_matrix */

/**
 * Print a square nxn matrix that is stored in vectorized format to
 * curses window.
 *
 * @param fp	FILE pointer
 * @param mat	vectorized matrix
 * @param n	dimension
 * @param row	row or column vectorization
 */
void wprint_vectorized_sq_matrix(WINDOW *wp, double *mat, size_t n, int row)
{
	size_t i, j;

	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if ((row ? mat[i*n + j] : mat[j*n + i]) < 1e-2)
				wprintw(wp, " %8.2e",
					row ? mat[i*n + j] : mat[j*n + i]);
			else
				wprintw(wp, " %8.3f",
					row ? mat[i*n + j] : mat[j*n + i]);
		}
		wprintw(wp, "\n");
	}
} /* wprint_vectorized_sq_matrix */

/**
 * Print vectorized matrix to file handle.
 *
 * @param fp	FILE handle
 * @param mat	vectorized matrix
 * @param n	number of rows
 * @param l	number of columns
 * @param row	row vectorized
 */
void fprint_vectorized_matrix(FILE *fp, double *mat, size_t n, size_t l,
	int row)
{
	size_t i, j;

	for (i = 0; i < n; ++i) {
		fprintf(fp, "%3lu", i);
		for (j = 0; j < l; ++j) {
			if ((row ? mat[i*l + j] : mat[j*l + i]) < 1e-2)
				fprintf(fp, " %8.2e",
					row ? mat[i*l + j] : mat[j*l + i]);
			else
				fprintf(fp, " %8.3f",
					row ? mat[i*l + j] : mat[j*l + i]);
		}
		fprintf(fp, "\n");
	}
} /* fprint_vectorized_matrix */

/**
 * Print matrix to file handle.
 *
 * @param fp	FILE handle
 * @param mat	matrix
 * @param n	number of rows
 * @param p	number of columns
 */
void fprint_matrix(FILE *fp, double **mat, size_t n, size_t p)
{
	size_t i, j;

	for (i = 0; i < n; ++i) {
		fprintf(fp, "%3lu", i);
		for (j = 0; j < p; ++j) {
			if (mat[i][j] < 1e-2)
				fprintf(fp, " %8.2e", mat[i][j]);
			else
				fprintf(fp, " %8.3f", mat[i][j]);
		}
		fprintf(fp, "\n");
	}
} /* fprint_matrix */

/**
 * Print vectorized matrix to curses window.
 *
 * @param fp	FILE handle
 * @param mat	vectorized matrix
 * @param n	number of rows
 * @param l	number of columns
 * @param row	row vectorized
 */
void wprint_vectorized_matrix(WINDOW *wp, double *mat, size_t n, size_t l,
	int row) {
	size_t i, j;

	for (i = 0; i < n; ++i) {
		wprintw(wp, "%3u", i);
		for (j = 0; j < l; ++j) {
			if ((row ? mat[i*l + j] : mat[j*l + i]) < 1e-2)
				wprintw(wp, " %8.2e",
					row ? mat[i*l + j] : mat[j*l + i]);
			else
				wprintw(wp, " %8.3f",
					row ? mat[i*l + j] : mat[j*l + i]);
		}
		wprintw(wp, "\n");
	}
} /* wprint_vectorized_matrix */

/**
 * Print vectorized matrix (unsigned int) to file handle.
 *
 * @param fp	FILE handle
 * @param mat	vectorized matrix
 * @param n	number of rows
 * @param l	number of columns
 * @param row	row vectorized
 */
void fprint_vectorized_uintmatrix(FILE *fp, unsigned int *mat, unsigned int n, unsigned int l,
	int row)
{
	size_t i, j;

	for (i = 0; i < n; ++i) {
		for (j = 0; j < l; ++j) {
			if ((row ? mat[i*l + j] : mat[j*l + i]) < 1e-2)
				fprintf(fp, " %u",
					row ? mat[i*l + j] : mat[j*l + i]);
			else
				fprintf(fp, " %u",
					row ? mat[i*l + j] : mat[j*l + i]);
		}
		fprintf(fp, "\n");
	}
} /* fprint_vectorized_uintmatrix */

int fread_size_ts(FILE *fp, size_t *idx, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		if (fscanf(fp, "%lu", &idx[i]) != 1)
			return(mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"expecting %u elements in index file\n", n));

	return NO_ERROR;
} /* fread_size_ts */

int fread_uints(FILE *fp, unsigned int *id, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		if (fscanf(fp, "%u", &id[i]) != 1)
			return(mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"expecting %u elements in index file\n", n));

	return NO_ERROR;
} /* fread_uints */

int read_uints(char const * const filename, unsigned int *id, size_t nreads)
{
	FILE *fp = fopen(filename, "r");

	if (!fp)
		return(mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename));

	int err = fread_uints(fp, id, nreads);

	fclose(fp);
	return err;
} /* read_uints */
