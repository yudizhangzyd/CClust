/**
 * @file cluster_lloyds.c
 * @author Yudi Zhang
 *
 * Functions for Lloyd's algorithm applied to fastq data assuming literal
 * quality scores and uniform substitution probabilties.  It also contains
 * the general interface to Lloyd's algorithm that will presumably one
 * day replace the code in kmodes.c::kmodes_lloyd().
 *
 */

#include "cluster_lloyds.h"
#include "khaplotype_data.h"


#include "array.h"
#include "error.h"
#include "math.h"
#include "kmodes.h"

/**
 * Find the number of unique categories in a data matrix arranged in single
 * vector of length n.  The data are (fundamentally) integers, and this
 * algorithm does not merely return max data.  It does not count unobserved
 * values.
 *
 * @param data		pointer to vectorized data
 * @param n		length of vectorized data
 * @param n_unique	number of unique categories (TBD)
 * @return		vector of unique category
 */
data_t *find_unique(data_t *data, unsigned int n, data_t *n_unique)
{
	*n_unique = 0;
	data_t *unique = NULL;
	unsigned int *count = NULL;

	CMAKE_1ARRAY(count, n);
	MAKE_1ARRAY(unique, n);

	for (unsigned int i = 0; i < n; ++i) {
		data_t ind = data[i];
		if (!count[ind])
			unique[(*n_unique)++] = ind;
		++count[ind];
	}

	data_t *ptr = realloc(unique, *n_unique);
	if (!ptr)
		return NULL;
	unique = ptr;
	qsort(unique, *n_unique, sizeof(*unique), compare);

	FREE_VECTOR(count);

	return unique;
} /* find_unique */

int compare(const void *a, const void *b)
{
	return ( *(data_t *)a - *(data_t *)b );
} /* compare */

/**
 * Lloyd's algorithm.  Alternate between updating cluster assignments and
 * cluster centers until convergence.
 *
 * @param d		void pointer data structure: defined by user
 * @param seeds		seeds
 * @param nclass	count in clusters (1 x k)
 * @param ic1		cluster assignments (1 x n)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param max_iter	maximum allowed iterations
 * @param cost		cost per cluster
 * @param ifault	pointer to error status
 * @param iter		pointer to number of iterations
 * @param step1		function pointer: 1st step of lloyds
 * @param step2		function pointer: 2nd step of lloyds
 * @return optimized criterion
 */
double
cluster_lloyds (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int *nclass,
	unsigned int *ic1,
	unsigned int n, unsigned int p, unsigned int K,
	unsigned int max_iter,
	double *cost,
	int *ifault,
	unsigned int *iter,
	lloyd_step1 step1, lloyd_step2 step2, compute_rule compute_sum_cost)
{
	unsigned int m = 0;
	int keep_going = 1;
	*ifault = KMODES_NO_ERROR;

	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KMODES_CALLER_INPUT_ERROR;
		/* [KSD, TODO] Let caller exit gracefully. */
		mmessage(ERROR_MSG, INVALID_USER_INPUT,
				 "Request too few or many cluasters");
		exit(EXIT_FAILURE); // return INFINITY;
	}

	/* initialize cluster counts */
	SETZERO_1ARRAY(nclass, K);

	/* assign to initial clusters and count things */
	step1(d, auxd, seeds, K, n, p, nclass, ic1, 1);

	step2(d, auxd, seeds, K, n, p, ic1);
	
	while (m++ < max_iter && keep_going) {
	
		keep_going = step1(d, auxd, seeds, K, n, p, nclass, ic1, 0);

		if (!keep_going)
			break;
		/* update haplotypes */
		step2(d, auxd, seeds, K, n, p, ic1);

	}
	
	if (m > max_iter)
		*ifault = KMODES_EXCEED_ITER_WARNING;

	*iter = m;

	return compute_sum_cost(d, auxd, seeds, ic1, cost, K, n, p);
}/* cluster_lloyds */

/**
 * Compute cost criterion. [PARALLEL]
 *
 * @param d		data structure
 * @param ic		cluster assignments (n x 1)
 * @param K		number of clusters
 * @param n		number of observations
 * @param p		number of coordinates
 * @return		return sum of cluster costs
 */

double compute_criterion_hap(void *d, void *auxd, data_t **seeds, unsigned int *ic, double *criterion, unsigned int K, unsigned int n, unsigned int p)
{
	data *dp = (data *)d;
	auxvar *aux_dp = (auxvar *) auxd;

	/* reset per-cluster cost */
	for (unsigned int k = 0; k < K; ++k)
		criterion[k] = 0.;

	/* compute per-cluster cost */
	for (unsigned int i = 0; i < n; ++i)	/* [PARALLEL] */
		for (unsigned int j = 0; j < p; ++j) {
			criterion[ic[i]] += dp->dmat[i][j] == seeds[ic[i]][j] ? aux_dp->prob_t[i][j] : aux_dp->prob_f[i][j];
		}

	/* compute total cost */
	double sum = 0;
	for (unsigned int k = 0; k < K; ++k)
		sum += criterion[k];

	return sum;
 
} /* compute_criterion */

/**
 * Function to carry out step 1 for fastq data using uniform substitution
 * probabilities and literal quality scores.  This step reassigns observations
 * to clusters.
 *
 * @param d		void pointer to the data structure
 * @param K		no. clusters
 * @param n		no. observations
 * @param p		no. coordinates
 * @param nclass	cluster counts (1 x K)
 * @param ic1		cluster assignments (1 x n)
 * @param do_init	is this an initialization step
 * @return keep_going
 */
int
fastq_lloyds_step1 (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic1,
	int do_init)
{

	unsigned int i, j, l, k;
	double maxll, ll;
	int keep_going = 0;
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;

	for (i = 0; i < n; i++) {
		/* identify closest haplotype */
		maxll = -INFINITY;
		k = -1;
		for (l = 0; l < K; ++l) {	/* compare to each haplotype */
			ll = 0;
			for (j = 0; j < p; ++j)
				ll += dp->dmat[i][j] == seeds[l][j]
					? aux_dp->prob_t[i][j] : aux_dp->prob_f[i][j];

			if (ll > maxll) {
				k = l;
				maxll = ll;
			}
		}
		if (do_init) {
			nclass[k]++;
			ic1[i] = k;
		} else if (ic1[i] != k) {
			keep_going = 1;
			nclass[ic1[i]]--;
			nclass[k]++;
			ic1[i] = k;
		}
	}
	
	return keep_going;
}/* fastq_lloyds_step1 */

/**
 * Step 2 for lloyd applied to fastq data using uniform substitution
 * probabilities and literal quality scores.  This step updates the
 * haplotype nucleotides.
 *
 * @param d		void pointer
 * @param K		no. cluster
 * @param n		no. observation
 * @param p		no. coordinate
 * @param ic1		cluster assignments (1 x n)
 */
void
fastq_lloyds_step2 (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *ic1)
{

	unsigned int i, j, k, l;
	double llj, max_llj;
	unsigned int len_index;
	data *dp = (data *)d;
	auxvar *aux_dp = (auxvar *) auxd;

	/* update haplotypes */
	for (k = 0; k < K; ++k) {

		/* store indices of data in cluster k */
		len_index = 0;
		for (i = 0; i < n; ++i)
			if (ic1[i] == k)
				aux_dp->index_in_cluster[len_index++] = i;

		for (j = 0; j < p; ++j) {
			// compute maxll of current haplotype
			max_llj = -INFINITY;
//			for (i = 0; i < len_index; ++i)
//				max_llj += dp->dmat[dp->index[i]][j]
//					== seeds[k][j]
//					? dp->prob_t[dp->index[i]][j]
//					: dp->prob_f[dp->index[i]][j];
			/* compute ll for all other possible nucleotides */
			for (l = 0; l < dp->tot_n_categories; ++l) {
//				if (dp->categories[l] == seeds[k][j])
//					continue;
				llj = 0;
				for (i = 0; i < len_index; ++i)
					llj += dp->dmat[aux_dp->index_in_cluster[i]][j]
						== dp->categories[l]
						? aux_dp->prob_t[aux_dp->index_in_cluster[i]][j]
						: aux_dp->prob_f[aux_dp->index_in_cluster[i]][j];
				
				if (llj > max_llj) {
					
					seeds[k][j] = dp->categories[l];
					max_llj = llj;
				}
			}
		}
	}
} /* fastq_lloyds_step2 */
