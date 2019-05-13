/**
 * Initialization file shared by kmodes and khaplotypes
 * @file init.c
 * @author Karin S. Dorman
 *
 */

#include "init.h"

static int compare_data_to_seed(data_t **x, unsigned int i, data_t *seed, unsigned int p);
/* k-modes initialization routines */
static int kmodes_init_random_seeds(data_t **x, unsigned int n, unsigned int p, unsigned int K, unsigned int k1, data_t **seeds, unsigned int *sd_idx);

/**
 * Randomly choose seeds from a known partition.  If options::K is larger than
 * the number of partitions, uniformly sample partitions, then distinct seeds
 * from within
 *
 * @param x		data (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of seeds to select
 * @param seeds		place to put seeds (K x p)
 * @param sd_idx	place to store selected seed indices (1 x K)
 * @param id		partition (1 x n)
 * @return		error status
 */
int
kmodes_init_random_from_partition (data_t **x, unsigned int n, unsigned int p,
				   unsigned int K, data_t **seeds, unsigned int *sd_idx, unsigned int *id)
{
	unsigned int true_k = 0;
	unsigned int *nc = NULL;
	unsigned int *sidx = NULL;
	int same;
	
	if (n < K)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Requesting %u clusters with only %u reads.\n", K, n);
	
	if (!sd_idx) {
		sidx = malloc(K * sizeof *sidx);
		
		if (!sidx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"seed index");
	} else
		sidx = sd_idx;
	
	/* determine the true number of clusters */
	for (unsigned int i = 0; i < n; ++i)
		if (id[i] > true_k)
			true_k = id[i];
	++true_k;
	
	/* count the observations in each cluster */
	nc = calloc(true_k, sizeof *nc);
	if (!nc)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "nc");
	
	for (unsigned int i = 0; i < n; ++i)
		++nc[id[i]];
	
	/* more than K clusters */
	if (K >= true_k) {
		
		/* pick one seed from each cluster */
		for (unsigned int k = 0; k < true_k; ++k) {
			unsigned int j = (unsigned int)((double) rand()
							/ RAND_MAX * nc[k]);
			sidx[k] = j;
			memcpy(seeds[k], x[j], p * sizeof **x);
		}
		
		/* remainding seeds */
		for (unsigned int k = true_k; k < K; ++k) {
			/* select cluster with replacement */
			unsigned int m = (unsigned int)((double) rand()
							/ RAND_MAX * K);
			if (m == K) --m;
			
			/* select distinct observation */
			do {
				same = 0;
				sidx[k] = (unsigned int) ((double) rand()
							  / RAND_MAX * nc[k]);
				if (sidx[k] == nc[k]) --sidx[k];
				for (unsigned l = 0; l < sidx[k]; ++l)
					if (sidx[l] == sidx[k] ||
					    !compare_data(x, sidx[l],
							  sidx[k], p)) {
						    same = 1;
						    break;
					    }
				memcpy(seeds[k], x[sidx[k]], p * sizeof **x);
			} while (same);
		}
	} else if (K < true_k) {
		/* select clusters without replacement */
		unsigned int *sel_k = malloc(K * sizeof *sel_k);
		sample(true_k, K, sel_k);
		for (unsigned int k = 0; k < K; ++k) {
			unsigned int j = (unsigned int)((double) rand()
							/ RAND_MAX * nc[sel_k[k]]);
			if (sidx) sidx[k] = j;
			memcpy(seeds[k], x[j], p * sizeof **x);
		}
		free(sel_k);
	}
	
	free(nc);
	
	return NO_ERROR;
} /* kmodes_init_random_from_partition */


/**
 * Initialize by randomly selecting seeds from a set.
 *
 * @param K		number of seeds to select
 * @param p		number of coordinates
 * @param n_ss		number of seeds in seedset
 * @param seeds		seeds to set (K x p)
 * @param seedset	choices to select from (n_ss x p)
 * @return		error status
 */
int
kmodes_init_random_from_set (unsigned int K, unsigned int p,
			     unsigned int n_ss, data_t **seeds, data_t **seedset)
{
	
	if (K >= n_ss)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Need seed set "
				"size (%u) to exceed K=%u\n", n_ss, K);
	
	/* if triggered check p*sizeof(data_t) can store uint & use pointers */
	if (n_ss > pow(2, 8*sizeof(data_t)))
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "");
	
	unsigned int k = 0, t = 0;
	
	while (k < K) {
		double u = rand() / (RAND_MAX + 1.);
		
		if ( (n_ss - t) * u >= K - k )
			++t;
		else	/* [TODO] assumes K fits in data_t */
			seeds[k++][0] = (data_t) t++;
	}
	
	for (k = 0; k < K; ++k)
		memcpy(seeds[k], seedset[(size_t) seeds[k][0]], p * sizeof **seeds);
	
	return NO_ERROR;
} /* kmodes_init_random_from_set */


/**
 * @param x		data (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param k1		number of seeds already selected
 * @param seeds		seeds (to be set)
 * @param sd_idx	indices of chosen observations
 * @return		error status
 */
static int
kmodes_init_random_seeds (data_t **x, unsigned int n, unsigned int p,
			  unsigned int K, unsigned int k1, data_t **seeds, unsigned int *sd_idx)
{
	unsigned int *sidx = NULL;
	int same;
	
	if (n < K)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Requesting %u clusters with only %u reads.\n", K, n);
	
	if (!sd_idx) {
		sidx = malloc(K * sizeof *sidx);
		
		if (!sidx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"seed index");
	} else
		sidx = sd_idx;
	for (unsigned j = 0; j < k1; ++j) {
		sidx[j] = n;
		for (unsigned int i = 0; i < n; ++i)
			if (!compare_data_to_seed(x, i, seeds[j], p)) {
				sidx[j] = i;
				break;
			}
	}
	
	for (unsigned int i = k1; i < K; ++i) {
		
		do {
			same = 0;
			sidx[i] = (unsigned int) ((double) rand() / RAND_MAX * n);
			for (unsigned int j = 0; j < i; ++j)
			/* sample without replacement */
			/* check for same or equal seeds */
				if (sidx[i] == sidx[j] || (sidx[j] < n &&
							   !compare_data(x, sidx[i], sidx[j], p))) {
					same = 1;
					break;
				}
		} while (same);
		memcpy(seeds[i], x[sidx[i]], p * sizeof **x);
	}
	
	if (!sd_idx)
		free(sidx);
	
	return NO_ERROR;
} /* kmodes_init_random_seeds */

/**
 * Returns human-friendly description of k-modes and k-haplotypes initializatiom method.
 *
 * @param method	method
 * @return		string constant
 */
const char *
kcluster_init_method (int method)
{
	if (method == KMODES_INIT_RANDOM_SEEDS)
		return "random seeds";
	else if (method == KMODES_INIT_RANDOM_FROM_PARTITION)
		return "random initialization from true partition";
	else if (method == KMODES_INIT_RANDOM_FROM_SET)
		return "random initialization from set of seeds";
	else
		return NULL;
}/* kcluster_init_method */

/**
 * Initialize k-modes and khaplotype.
 *
 * @param x		data (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param k		number of clusters
 * @param k1		number of seeds deterministically selected
 * @param seeds		initial modes to be chosen by this function
 * @param sidx		seed indices to be set (if initializer's seeds are observations)
 * @param method 	method of initialization (see kmodes.h)
 * @return		error status
 */
int
kcluster_init (data_t **x, unsigned int n, unsigned int p, unsigned int k, unsigned int k1, data_t **seeds, unsigned int *sidx, int method)
{
	if (method == KMODES_INIT_RANDOM_SEEDS) {
		return kmodes_init_random_seeds(x, n, p, k, k1, seeds, sidx);
	} else
		return mmessage(ERROR_MSG, KMODES_INVALID_INITIALIZATION_METHOD,
				kmodes_error(KMODES_INVALID_INITIALIZATION_METHOD));
	return NO_ERROR;
} /* kmodes_init */

static int
compare_data_to_seed (
		      data_t **x,
		      unsigned int i,
		      data_t *seed,
		      unsigned int p)
{
	for (unsigned int j = 0; j < p; ++j) {
		if (x[i][j] > seed[j])
			return 1;
		else if (x[i][j] < seed[j])
			return -1;
	}
	return 0;
} /* compare_data_to_seed */

/**
 * Compare two observations to allow ordering and testing for equality.
 *
 * @param x	data (? x n)
 * @param i	index of first observation
 * @param j	index of second observation
 * @param n	number of coordinates
 * @return	-1, 0, 1 order indicator
 */
int
compare_data (
	      data_t **x,
	      unsigned int i,
	      unsigned int j,
	      unsigned int n)
{
	for (unsigned int l = 0; l < n; ++l) {
		if (x[i][l] > x[j][l])
			return 1;
		else if (x[i][l] < x[j][l])
			return -1;
	}
	return 0;
} /* compare_data */
