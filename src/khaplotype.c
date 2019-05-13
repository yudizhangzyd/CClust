/**
 * @file khaplotype.c
 * @author Yudi Zhang
 *
 * Functions for efficient algorithm applied to fastq data assuming literal
 * quality scores and uniform substitution probabilties.  It also contains
 * the general interface to Lloyd's algorithm that will presumably one
 * day replace the code in kmodes.c::kmodes_lloyd().
 *
 * [TODO] data::tot_n_categories: forces all sites to use the same categories
 * [TODO] make algorithm more flexible
 * - separate data structure elements: auxiliary variables needed for algorithm vs. data storage
 * - the algorithm knows too much about the data structure and the data knows too much about the algorithm
 * - separate them so that a person can run the algorithm on their own data structure without having to include elements for the auxiliary variables
 *
 */

#include "khaplotype.h"
#include "khaplotype_data.h"
#include "khaplotype_model.h"
#include "khaplotype_option.h"
#include "khaplotype_initialize.h"
#include "array.h"

/* DEBUGGING: check where most updates are happening
unsigned int *update_j = NULL;
unsigned int *update_i = NULL;
*/

static void init_auxiliary_variables(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic1);
static void init_auxiliary_variables_hw(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic1);
static void init_haplotypes(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p);
static void move_observation(void *d, void *auxd, data_t **seeds, unsigned int *nclass, unsigned int n, unsigned int i, unsigned int from, unsigned int to);
//static void init_haplotypes(void *d, data_t **seeds, unsigned int K, unsigned int n, unsigned int p);

/**
 * Compute cost criterion.
 *
 * @param d           data structure
 * @param ic          cluster assignments (n x 1)
 * @param criterion   cost per cluster
 * @param K           number of clusters
 * @param n           number of observations
 * @return            return sum of cluster costs
 */
double
compute_cost (
	void *d,
	unsigned int *ic,
	double *criterion,
	unsigned int K,
	unsigned int n)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
#endif
	double sum = 0;
	auxvar *dp = (auxvar *) d;

	SETZERO_1ARRAY(criterion, K);

#ifdef DEBUGGING_CODE
	debug_msg(DEBUG_I, fxn_debug, "ic: ");
	if (fxn_debug >= DEBUG_I)
		fprint_uints(stderr, ic, n, 1, 1);
#endif
	
	for (unsigned int i = 0; i < n; ++i) {
#ifdef DEBUGGING_CODE
		debug_msg(DEBUG_I, fxn_debug, "v[%u]: ", i);
		if (fxn_debug >= DEBUG_I)
			fprint_doubles(stderr, dp->v_ik[i], K, 3, 1);
#endif
		criterion[ic[i]] += dp->v_ik[i][ic[i]];
	}
	
	for (unsigned int k = 0; k < K; ++k)
		sum += criterion[k];
#ifdef DEBUGGING_CODE
	debug_msg(DEBUG_I, fxn_debug, "Log likelihood in each cluster.");
	fprint_doubles(stderr, criterion, K, 3, 0);
	fprintf(stderr, ": %f\n", sum);
#endif
	
	return sum;
} /* compute_cost */

/**
 * Step 2 for lloyd applied to fastq data using uniform substitution
 * probabilities and literal quality scores.  This step updates the
 * haplotype nucleotides.
 *
 * @param d		void pointer
 * @param seeds		seeds
 * @param K		no. cluster
 * @param n		no. observation
 * @param p		no. coordinate
 * @param ic1		cluster assignments (1 x n)
 * @param do_init	if it is at the initial phase
 */
void
fastq_lloyds_efficient_step2 (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *ic1,
	int do_init)
{
	
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int i, k, j, l, len_index;
	double e_kjN_diff;
	
	/* Update e_kjn */
	if (!do_init)
		for (i = 0; i < n; ++i) {
			if (aux_dp->last_assign[i] == ic1[i])
				continue;
			for (j = 0; j < p; ++j)
				for (l = 0; l < dp->tot_n_categories; ++l) {
					e_kjN_diff =
						dp->dmat[i][j]
						== dp->categories[l]
						? aux_dp->prob_t[i][j]
						: aux_dp->prob_f[i][j];
					aux_dp->e_kjN[aux_dp->last_assign[i]][j][l]
						-= e_kjN_diff;
					aux_dp->e_kjN[ic1[i]][j][l] += e_kjN_diff;
				}
		}
	else
		for (k = 0; k < K; ++k)
			for (j = 0; j < p; ++j)
				for (l = 0; l < dp->tot_n_categories; ++l)
					aux_dp->e_kjN[k][j][l] = 0;
	
	/* Store e_kjN; Second step: recompute the center */
	for (k = 0; k < K; ++k) {
		len_index = 0;
		
		if (do_init)
			for (i = 0; i < n; i++)
				if (ic1[i] == k)
					aux_dp->index_in_cluster[len_index++] = i;
		
		for (j = 0; j < p; ++j) {
			/* Store the current center for the following comparision */
			aux_dp->last_cent[k][j] = seeds[k][j];
			
			double max_llj = -INFINITY;
			int max_id = -1;
			
			for (l = 0; l < dp->tot_n_categories; ++l) {
				if (do_init)
					for (i = 0; i < len_index; ++i)
						aux_dp->e_kjN[k][j][l] +=
							dp->dmat[aux_dp->index_in_cluster[i]][j]
							== dp->categories[l]
							? aux_dp->prob_t[aux_dp->index_in_cluster[i]][j]
							: aux_dp->prob_f[aux_dp->index_in_cluster[i]][j];
				/* find the largest e_kjN[k][j][l] for each site */
				if (aux_dp->e_kjN[k][j][l] > max_llj) {
					max_llj = aux_dp->e_kjN[k][j][l];
					max_id = l;
				}
				
			}
			seeds[k][j] = dp->categories[max_id]; //Compute new Hk
		}
	}
}/* fastq_lloyds_efficient_step2 */

/**
 * Function to carry out step 1 for fastq data using uniform substitution
 * probabilities and literal quality scores.  This step reassigns observations
 * to clusters.
 *
 * @param d		void pointer to the data structure
 * @param seeds		seeds
 * @param K		no. clusters
 * @param n		no. observations
 * @param p		no. coordinates
 * @param nclass	cluster counts (1 x K)
 * @param ic1		cluster assignments (1 x n)
 * @param do_init	is this an initialization step
 * @return		keep_going
 */
int
fastq_lloyds_efficient_step1 (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic1,
	int do_init)
{
	
	unsigned int i, k, j;
	int keep_going = 0;
	double v_ik_diff;
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;

	if (do_init)
		for (i = 0; i < n; ++i)
			for (k = 0; k < K; ++k)
				aux_dp->v_ik[i][k] = 0;
	
	/* Reassign reads to new centers, detect change in centers */
	for (i = 0; i < n; ++i) {        /* for each observation */

		/* Store the last cluster assignment */
		if (!do_init)
			aux_dp->last_assign[i] = ic1[i];

		double max = -INFINITY;
		unsigned int max_id = -1;
		unsigned int last_assign = aux_dp->last_assign[i];
		
		for (k = 0; k < K; ++k) {    /* compare to each haplotype */
			for (j = 0; j < p; ++j) {
				if (do_init) {
					aux_dp->v_ik[i][k] +=
						dp->dmat[i][j] == seeds[k][j]
						? aux_dp->prob_t[i][j]
						: aux_dp->prob_f[i][j];
				
				} else if (aux_dp->last_cent[k][j] != seeds[k][j]) {
					
					v_ik_diff = aux_dp->prob_t[i][j]
							- aux_dp->prob_f[i][j];
					
					if (dp->dmat[i][j] == seeds[k][j])
						aux_dp->v_ik[i][k] += v_ik_diff;
					else if (dp->dmat[i][j]
						== aux_dp->last_cent[k][j])
						aux_dp->v_ik[i][k] -= v_ik_diff;
				}
			}
			if (aux_dp->v_ik[i][k] > max) {
				max = aux_dp->v_ik[i][k];
				max_id = k;
			}
		}

		if (do_init) {
			nclass[max_id]++;
			ic1[i] = max_id;
		} else if (last_assign != max_id) {
			keep_going = 1;
			nclass[max_id]++;
			nclass[last_assign]--;
			ic1[i] = max_id;
		}
	}

	return keep_going;
} /* fastq_lloyds_efficient_step1 */


/**
 * Efficient version of Lloyd's algorithm.  Alternate between updating cluster
 * assignments and cluster centers until convergence.
 *
 * @param d void pointer data structure: defined by user
 * @param seeds seeds
 * @param nclass count in clsuters (1 x k)
 * @param ic1 cluster assignments (1 x n)
 * @param n number of observations
 * @param p number of coordinates
 * @param K number of clusters
 * @param max_iter maximum allowed iterations
 * @param cost cost per cluster
 * @param ifault pointer to error status
 * @param iter pointer to number of iterations
 * @param step1 function pointer: user defined way to carry out the 1st step of lloyds
 * @param step2 function pointer: user defined way to carry out the 2nd step of lloyds
 * @return optimized criterion
 */
double
cluster_lloyds_efficient (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int *nclass,
	unsigned int *ic1,
	unsigned int n, unsigned int p, unsigned int K,
	unsigned int max_iter,
	double *cost,
	int *ifault, unsigned int *iter,
	lloyd_step1 step1,
	modified_lloyd_step2 step2,
	compute_crit compute_sum_cost)
{
	
	*ifault = KHAPLOTYPE_NO_ERROR;
	unsigned int m = 0;
	int keep_going = 1;
	
	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KHAPLOTYPE_CALLER_INPUT_ERROR;
		mmessage(ERROR_MSG, INVALID_USER_INPUT,
			 "Request too few or many cluasters");
		exit(EXIT_FAILURE); // return INFINITY;
	}
	
	/* initialize cluster counts and cost */
	SETZERO_1ARRAY(nclass, K);

	step1(d, auxd, seeds, K, n, p, nclass, ic1, 1);
	step2(d, auxd, seeds, K, n, p, ic1, 1);
	
	while (m++ < max_iter && keep_going) {
		
		keep_going = step1(d, auxd, seeds, K, n, p, nclass, ic1, 0);
		
		if (!keep_going)
			break;
		
		step2(d, auxd, seeds, K, n, p, ic1, 0);
	}
	
	if (m > max_iter)
		*ifault = KHAPLOTYPE_EXCEED_ITER_WARNING;
	
	*iter = m;
	
	return compute_sum_cost(auxd, ic1, cost, K, n);
	
} /* cluster_lloyds_efficient */

/**
 * Modified MacQueen's algorithm.  Initialization step is the same as Lloyd's
 * efficient.
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
 * @return		optimized criterion
 */
double
cluster_macqueen (
	void *d,
	void *auxd,
	unsigned int n, unsigned int p,
	data_t **seeds,
	unsigned int K,
	unsigned int *ic1,
	unsigned int *nclass,
	unsigned int max_iter,
	double *cost,
	int *ifault,
	unsigned int *iter,
	macq_ini ini,
	macq_iter itr,
	compute_crit compute_sum_cost)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
#endif
	*ifault = KHAPLOTYPE_NO_ERROR;
	unsigned int k, m = 0;
	int keep_going = 1;

#ifdef DEBUGGING_CODE
if (fxn_debug >= DEBUG_I) {
	debug_msg(DEBUG_I, fxn_debug, "Haplotypes:\n");
	for (unsigned int i = 0; i < K; ++i) {
		for (unsigned int j = 0; j < p; ++j)
			fprintf(stderr, " %u", seeds[i][j]);
		fprintf(stderr, "\n");
	}
}
#endif
	
	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KHAPLOTYPE_CALLER_INPUT_ERROR;
		mmessage(ERROR_MSG, INVALID_USER_INPUT,
			 "Request too few or many clusters");
		exit(EXIT_FAILURE); // return INFINITY;
	}
	
	ini(d, auxd, seeds, K, n, p, nclass, ic1);
	
	/* Iteration */
	while (m++ < max_iter && keep_going) {
		keep_going = itr(d, auxd, seeds, K, n, p, nclass, ic1);
		if (!keep_going)
			break;
	}

	m--;
	*iter = m;
	
	if (m > max_iter && max_iter)
		*ifault = KHAPLOTYPE_EXCEED_ITER_WARNING;
	
	/* error if initialization routine produces empty cluster */
	for (k = 0; k < K; ++k)
		if (nclass[k] == 0) {
			*ifault = KHAPLOTYPE_NULL_CLUSTER_ERROR;
			return INFINITY;
		}
	
	return compute_sum_cost(auxd, ic1, cost, K, n);
} /* cluster_macqueen */


/**
 * Iterate fastq Macqueen.
 *
 * @param d		void pointer to data object
 * @param seeds		current centers
 * @param K		number of clusters
 * @param n		number of observations
 * @param p		dimension of observations
 * @param nclass	size of clusters
 * @param ic1		cluster assignments of observations
 * @return		whether to keep iterating
 */
int
fastq_macqueen_iter (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic1)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
#endif
	unsigned int i, k, j, l;
	unsigned int lassign;
	int keep_going = 0;
	double v_ik_diff, e_kjN_diff;
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	
	for (i = 0; i < n; ++i) {
		double max = -INFINITY;
		unsigned int max_id = 0;

		aux_dp->last_assign[i] = ic1[i];
		lassign = ic1[i];
		
		/* Update cluster assignment */
		for (k = 0; k < K; ++k)
			if (aux_dp->v_ik[i][k] > max) {
				max = aux_dp->v_ik[i][k];
				max_id = k;
			}

		if (lassign == max_id)
			continue;

		/* cluster assignment changed */
		keep_going = 1;
#ifdef DEBUGGING_CODE
		debug_msg(DEBUG_I, fxn_debug, "moving %u from %u to %u\n",
							i, lassign, max_id);
#endif
		/* update cluster counts */
		nclass[max_id]++;
		nclass[lassign]--;
		ic1[i] = max_id;

		/* [KSD] Could speed up to only update changed sites. */
		/* [YZ] Isn't it only updating the sites? Implemented at "dp->last_cent[lassign][j] != seeds[lassign][j]". */
		for (j = 0; j < p; ++j) {
			
			double max_lst_llj = -INFINITY, max_cur_llj = -INFINITY;
			int max_lst_id = -1, max_cur_id = -1;

			for (l = 0; l < dp->tot_n_categories; ++l) {

				e_kjN_diff =
					dp->dmat[i][j] == dp->categories[l]
					? aux_dp->prob_t[i][j] : aux_dp->prob_f[i][j];
				aux_dp->e_kjN[lassign][j][l] -= e_kjN_diff;
				aux_dp->e_kjN[max_id][j][l] += e_kjN_diff;

				/* find largest e_kjN[k][j][l] for each site */
				if (aux_dp->e_kjN[lassign][j][l] > max_lst_llj) {
					max_lst_llj = aux_dp->e_kjN[lassign][j][l];
					max_lst_id = l;
				}

				if (aux_dp->e_kjN[max_id][j][l] > max_cur_llj) {
					max_cur_llj = aux_dp->e_kjN[max_id][j][l];
					max_cur_id = l;
				}
			}
			/* Update haplotypes */
			seeds[lassign][j] = dp->categories[max_lst_id];
			seeds[max_id][j] = dp->categories[max_cur_id];
			
			/* Update vik. Last center should be the one used to
			 * compute the last vik, so after iterate i from 1 to n,
			 * we update the values of last_cent
			 */
			if (aux_dp->last_cent[lassign][j] != seeds[lassign][j]) {
#ifdef DEBUGGING_CODE
				debug_msg(DEBUG_I, fxn_debug, "update cluster "
					"%u haplotype at site %u\n", lassign, j);
#endif
				for (unsigned int i1 = 0; i1 < n; ++i1) {
					v_ik_diff = aux_dp->prob_t[i1][j]
							- aux_dp->prob_f[i1][j];
					if (dp->dmat[i1][j] == seeds[lassign][j])
						aux_dp->v_ik[i1][lassign]
								+= v_ik_diff;
					else if (dp->dmat[i1][j]
						== aux_dp->last_cent[lassign][j])
						aux_dp->v_ik[i1][lassign]
							-= v_ik_diff;
				}
				aux_dp->last_cent[lassign][j] = seeds[lassign][j];
			}

			if (aux_dp->last_cent[max_id][j] != seeds[max_id][j]) {
#ifdef DEBUGGING_CODE
				debug_msg(DEBUG_I, fxn_debug, "update cluster "
					"%u haplotype at site %u\n", max_id, j);
#endif
				for (unsigned int i1 = 0; i1 < n; ++i1) {
					v_ik_diff = aux_dp->prob_t[i1][j]
							- aux_dp->prob_f[i1][j];
					if (dp->dmat[i1][j] == seeds[max_id][j])
						aux_dp->v_ik[i1][max_id]
								+= v_ik_diff;
					else if (dp->dmat[i1][j]
						== aux_dp->last_cent[max_id][j])
						aux_dp->v_ik[i1][max_id]
								-= v_ik_diff;
				}
				aux_dp->last_cent[max_id][j] = seeds[max_id][j];
			}
		}
	}
	return keep_going;
} /* fastq_macqueen_iter */


/*
 * Initialize MacQueen.
 *
 * @param d		void pointer to data object
 * @param seeds		initial seeds
 * @param K		number of clusters
 * @param n		number of observations
 * @param p		number of coordinates
 * @param nclass	number of members in each cluster
 * @param ic1		cluster assignments
 */
void
fastq_macqueen_ini (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic1)
{

/* now done in init_auxiliary_variables()
	SETZERO_1ARRAY(nclass, K);
*/

	/* initialize auxiliary variables: ekjn, vik */
	init_auxiliary_variables(d, auxd, seeds, K, n, p, nclass, ic1);

	/* update haplotypes */
	init_haplotypes(d, auxd, seeds, K, n, p);

/* Now done in init_haplotypes:
	for (k = 0; k < K; ++k)
		for (j = 0; j < p; ++j)
			dp->last_cent[k][j] = seeds[k][j];
*/
} /* fastq_macqueen_ini */


/**
 * Hartigan and Wong's algorithm for clustering categorical data.
 *
 * @param d		void pointer data structure: defined by user
 * @param seeds		seeds
 * @param nclass	count in clsuters (1 x k)
 * @param ic		cluster assignments (1 x n)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param max_iter	maximum allowed iterations
 * @param cost		cost per cluster
 * @param ifault	pointer to error status
 * @param niter		pointer to number of iterations
 * @param init		function pointer: initialization
 * @param iter		function pointer: iterator
 * @returnr		optimized criterion
 **/
double
cluster_hw (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int *nclass,
	unsigned int *ic,
	unsigned int n, unsigned int p, unsigned int K,
	unsigned int max_iter,
	double *cost,
	int *ifault,
	unsigned int *niter,
	hw_init init,
	hw_iter iter)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
#endif
	
	*ifault = KHAPLOTYPE_NO_ERROR;
	unsigned int i, k, len_index, indx = 0;
//	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	
	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KHAPLOTYPE_CALLER_INPUT_ERROR;
		mmessage(ERROR_MSG, INVALID_USER_INPUT,
			 "Request too few or many cluasters");
		exit(EXIT_FAILURE); // return INFINITY;
	}

	*niter = 0;

	/* Initialization */
	init(d, auxd, seeds, K, n, p, nclass, ic, cost);
	
#ifdef DEBUGGING_CODE
	printf("cluster size:\n");
	fprint_uints(stdout, nclass, K, 4, 1);
	printf("assignment:\n");
	fprint_uints(stdout, ic, n, 2, 1);
#endif

	/* Optimal transfer stage */
	for (unsigned int m = 0; m < max_iter; m++) {
		++(*niter);
#ifdef DEBUGGING_CODE
		printf("iter %d\n", m);
#endif
		for (i = 0; i < n; ++i) {
			indx++;
			/* If in the live cluster, consider move to all the clusters */
			if (i < aux_dp->live[ic[i]]) { /* [KSD, TODO] We need to cleanup these unsigned int vs. int comparisons: one way is to drop to int types for n, p, and K OR change -= n below. */
				aux_dp->index = aux_dp->all_index;
#ifdef DEBUGGING_CODE
	printf("%d is in live set %d\n", i, ic[i]);
#endif
				iter(d, auxd, seeds, K, n, p, nclass, ic, cost, i, 1, &indx);
			}
			/* If not in the live set, only consider move to live sets */
			else {
				aux_dp->index = aux_dp->some_index;
				len_index = 0;
				/* Find the live set */
				for (k = 0; k < K; ++k)
					if (i < aux_dp->live[k])
						aux_dp->some_index[len_index++] = k;
					else if (k == ic[i])
						aux_dp->some_index[len_index++] = k;
#ifdef DEBUGGING_CODE
	printf("%d is in dead set %d\n", i, ic[i]);
	fprint_ints(stdout, aux_dp->live, len_index, 1, 1);
					
#endif

				iter(d, auxd, seeds, len_index, n, p, nclass, ic, cost, i, 0, &indx);
			}
			/* no transfers: algorithm has converged */
			if (indx == n)
				break;
		}
		
		for (k = 0; k < K; ++k)
			aux_dp->live[k] -= n;
		
		if (indx == n)
			break;
	}

	COPY_1ARRAY(cost, aux_dp->last_cost, K);
	double sum = 0;
	for (unsigned int k = 0; k < K; ++k)
		sum += cost[k];
	
	return sum;
}/* cluster_hw */


/**
 * Initialization for HW.
 *
 * @param d		void pointer to data structure: defined by user
 * @param seeds		seeds
 * @param nclass	count in clusters (1 x k)
 * @param ic		cluster assignments (1 x n)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param cost		cost per cluster
 */
void
hw_init_default (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic,
	double *cost)
{
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int i, k, j, l;
	
	/* initialize auxiliary variables: ekjn and vik */
	init_auxiliary_variables(d, auxd, seeds, K, n, p, nclass, ic);
	
	/* second step: update center sequences */
	init_haplotypes(d, auxd, seeds, K, n, p);
	
	/* initialize the cost */
	SETZERO_1ARRAY(cost, K);
	for (i = 0; i < n; ++i)
		cost[ic[i]] += aux_dp->v_ik[i][ic[i]];
	
	/* store last v_ik for easy reset */
	for (i = 0; i < n; ++i)
		for (k = 0; k < K; ++k)
			aux_dp->last_vik[i][k] = aux_dp->v_ik[i][k];

	for (k = 0; k < K; ++k) {
		aux_dp->live[k] = n + 1;	/* all clusters start to live */

		/* store the last cost */
		aux_dp->last_cost[k] = cost[k];

		for (j = 0; j < p; ++j) {
			/* Store the current center for the following comparision */
//			dp->last_cent[k][j] = seeds[k][j];

			/* store last ekjn for easy reset */
			for (l = 0; l < dp->tot_n_categories; ++l)
				aux_dp->last_ekjn[k][j][l] = aux_dp->e_kjN[k][j][l];
		}
	}
} /* hw_init_default */

/**
 * Optimal transfer stage

 * @param d		data pointer
 * @param seeds		centers
 * @param K		no. clusters to compare
 * @param n		no. observations
 * @param p		no. of coordinates
 * @param nclass	no. of observations in each cluster
 * @param ic		observation assignment
 * @param cost		cost
 * @param i		ith observation
 * @param is_live	to iniciate if ith observation is in the live set or not
 * @param indx		first 0, number of consecutive observations not transferred
 */
void
hw_iter_default (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic,
	double *cost,
	unsigned int i,
	int is_live,
	unsigned int *indx)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
#endif
	
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int i1, k, j, l;
	double max_llj, v_ik_diff;
	int max_id;
	double *diff = NULL;
	
	/* Store the last assignment */
	aux_dp->last_assign[i] = ic[i];
	
#ifdef DEBUGGING_CODE
	printf("old seeds:\n");
	fprint_matrix(stderr ,aux_dp->last_cent, K, p);
	printf("old v_ik:\n");
	fprint_matrix(stderr, aux_dp->last_vik, n, K);
#endif
	
	for (k = 0; k < K; ++k)
		/* move i to another cluster */
		if (aux_dp->index[k] != aux_dp->last_assign[i])
			for (j = 0; j < p; ++j)
				for (l = 0; l < dp->tot_n_categories; ++l) {
					if (dp->dmat[i][j] == dp->categories[l])
						aux_dp->e_kjN[aux_dp->index[k]][j][l]
							+= aux_dp->prob_t[i][j];
					else
						aux_dp->e_kjN[aux_dp->index[k]][j][l]
							+= aux_dp->prob_f[i][j];
				}
		/* move i from its current cluster */
		else
			for (j = 0; j < p; ++j)
				for (l = 0; l < dp->tot_n_categories; ++l) {
					if (dp->dmat[i][j] == dp->categories[l])
						aux_dp->e_kjN[aux_dp->index[k]][j][l]
							-= aux_dp->prob_t[i][j];
					else
						aux_dp->e_kjN[aux_dp->index[k]][j][l]
							-= aux_dp->prob_f[i][j];
				}

	/* Update the haplotypes and vik, to do: make it to take the function */
	for (k = 0; k < K; ++k)
		for (j = 0; j < p; ++j) {
			
			max_llj = -INFINITY;
			max_id = -1;
			
			for (l = 0; l < dp->tot_n_categories; ++l)
				if (aux_dp->e_kjN[aux_dp->index[k]][j][l] > max_llj) {  /*find the largest e_kjN[k][j][l] for each site*/
					max_llj = aux_dp->e_kjN[aux_dp->index[k]][j][l];
					max_id = l;
				}

			seeds[aux_dp->index[k]][j] = dp->categories[max_id]; //Compute new Hk
			
			if (seeds[aux_dp->index[k]][j] != aux_dp->last_cent[aux_dp->index[k]][j])
				
				for (i1 = 0; i1 < n; ++i1) {
					v_ik_diff = aux_dp->prob_t[i1][j] - aux_dp->prob_f[i1][j];
					if (dp->dmat[i1][j] == seeds[aux_dp->index[k]][j])
						aux_dp->v_ik[i1][aux_dp->index[k]] += v_ik_diff;
					else if (dp->dmat[i1][j] == aux_dp->last_cent[aux_dp->index[k]][j])
						aux_dp->v_ik[i1][aux_dp->index[k]] -= v_ik_diff;
				}
		}
	
#ifdef DEBUGGING_CODE
	printf("new seeds:\n");
	fprint_matrix(stdout, seeds, K, p);
	printf("new v_ik:\n");
	fprint_matrix(stdout, aux_dp->v_ik, n, K);
#endif
	
	/* recompute the cost, first initialize it */
	for (k = 0; k < K; ++k)
		cost[aux_dp->index[k]] = aux_dp->v_ik[i][aux_dp->index[k]];
	cost[aux_dp->last_assign[i]] = 0;
	
	if (is_live) {
		for (i1 = 0; i1 < n; ++i1)
			cost[ic[i1]] += aux_dp->v_ik[i1][ic[i1]];
		cost[ic[i]] -= aux_dp->v_ik[i][ic[i]];
	}
	else {
		for (k = 0; k < K; ++k)
			for (i1 = 0; i1 < n; ++i1)
				if (i1 != i && ic[i1] == aux_dp->index[k])
					cost[aux_dp->index[k]] += aux_dp->v_ik[i1][aux_dp->index[k]];
	}
#ifdef DEBUGGING_CODE
	printf("last cost:\n");
	fprint_doubles(stdout, aux_dp->last_cost, K, 3, 1);
	printf("cost:\n");
	fprint_doubles(stdout, cost, K, 3, 1);
#endif
	
	double max_dif = -INFINITY;
	MAKE_1ARRAY(diff, K);
	/* find the cluster that increases likelihood the most */
	double delta = 1e-6;
	double temp1 = 0.;
	double temp2 = 0.;
	double max = 0.;
	for (k = 0; k < K; ++k) {
		if (aux_dp->index[k] != aux_dp->last_assign[i]) {
			temp1 = cost[aux_dp->index[k]] + cost[aux_dp->last_assign[i]];
			temp2 = aux_dp->last_cost[aux_dp->index[k]]
				+ aux_dp->last_cost[aux_dp->last_assign[i]];
			diff[k] = temp1 - temp2;
			max = fabs(temp1);
			if (fabs(temp2) > fabs(temp1))
				max = fabs(temp2);
			if (fabs(diff[k]) <= (delta * max))
				diff[k] = 0.;
		} else {
			diff[k] = 0.;
		}
		if (diff[k] > max_dif) {
			max_dif = diff[k];
			ic[i] = aux_dp->index[k];
		}
	}
	free(diff);
	
	if (aux_dp->last_assign[i] != ic[i]) {
		*indx = 0;
		nclass[ic[i]]++;
		nclass[aux_dp->last_assign[i]]--;
		
		for (k = 0; k < K; ++k) {
			/* update the unaffected sets */
			if (aux_dp->index[k] != ic[i] && aux_dp->index[k] != aux_dp->last_assign[i]) {
				/* Reset vik and ekjn to the last value for the clusters that are not involved in the transformation of observation i */
				for (i1 = 0; i1 < n; ++i1)
					aux_dp->v_ik[i1][aux_dp->index[k]] = aux_dp->last_vik[i1][aux_dp->index[k]];
				for (j = 0; j < p; ++j) {
					for (l = 0; l < dp->tot_n_categories; ++l)
						aux_dp->e_kjN[aux_dp->index[k]][j][l] = aux_dp->last_ekjn[aux_dp->index[k]][j][l];
				}
			}
			/* update the affected sets */
			else
			{
				aux_dp->last_cost[aux_dp->index[k]] = cost[aux_dp->index[k]];
				aux_dp->live[aux_dp->index[k]] = n + i + 1;
				
				for (i1 = 0; i1 < n; ++i1)
				/* Store the last v_ik of all the ovservations for reset */
					aux_dp->last_vik[i1][aux_dp->index[k]] = aux_dp->v_ik[i1][aux_dp->index[k]];
				
				for (j = 0; j < p; ++j) {
					/* Store the current center for the following comparision */
					aux_dp->last_cent[aux_dp->index[k]][j] = seeds[aux_dp->index[k]][j];
					/* Store the last ekjn for reset */
					for (l = 0; l < dp->tot_n_categories; ++l)
						aux_dp->last_ekjn[aux_dp->index[k]][j][l] = aux_dp->e_kjN[aux_dp->index[k]][j][l];
				}
			}
		}
		
#ifdef DEBUGGING_CODE
		printf("move %d from %d to %d\n", i, aux_dp->last_assign[i], ic[i]);
		printf("new last_cost for the next round:\n");
		fprint_doubles(stdout, aux_dp->last_cost, K, 3, 1);
		printf("new v_ik for the next round:\n");
		fprint_matrix(stdout, aux_dp->v_ik, n, K);
#endif
		
	}
	else {
		/* Reset vik, ekjn and the centers to their last values */
		for (k = 0; k < K; ++k) {
			for (i1 = 0; i1 < n; ++i1)
				aux_dp->v_ik[i1][aux_dp->index[k]]
					= aux_dp->last_vik[i1][aux_dp->index[k]];
			for (j = 0; j < p; ++j)
				for (l = 0; l < dp->tot_n_categories; ++l)
					aux_dp->e_kjN[aux_dp->index[k]][j][l]
						= aux_dp->last_ekjn[aux_dp->index[k]][j][l];
		}
#ifdef DEBUGGING_CODE
		printf("Didn't move! Everything should be the same as considering moving %d\n", i);
		printf("new v_ik for the next round:\n");
		fprint_matrix(stdout, aux_dp->v_ik, n, K);
#endif
	}
} /* hw_iter_default */

/**
 * Hartigan and Wong's algorithm for clustering haplotypes (faster).
 *
 * @param d		void pointer data structure: defined by user
 * @param seeds		seeds
 * @param nclass	count in clusters (1 x K)
 * @param ic		cluster assignments (1 x n)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param max_iter	maximum allowed iterations
 * @param cost		cost per cluster
 * @param ifault	pointer to error status
 * @param niter		pointer to number of iterations
 * @param init		function pointer: initialization
 * @param iter		function pointer: iterator
 * @return optimized criterion
 **/
double
cluster_hw_fast (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int *nclass,
	unsigned int *ic,
	unsigned int n, unsigned int p, unsigned int K,
	unsigned int max_iter,
	double *cost,
	int *ifault,
	unsigned int *niter,
	hw_fast_init init,
	hw_fast_iter iter)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
#endif
	
	*ifault = KHAPLOTYPE_NO_ERROR;
	unsigned int i, j, k, len_index, indx = 0;
	double sum = 0;
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;

	*niter = 0;
	
	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KHAPLOTYPE_CALLER_INPUT_ERROR;
		mmessage(ERROR_MSG, INVALID_USER_INPUT,
			 "Request too few or many cluasters");
		exit(EXIT_FAILURE); // return INFINITY;
	}
	
	/* initialize cluster counts and cost */
	SETZERO_1ARRAY(cost, K);
	
	/* Initialization */
	init(d, auxd, seeds, K, n, p, nclass, ic);
	
#ifdef DEBUGGING_CODE
	printf("cluster size:\n");
	fprint_uints(stdout, nclass, K, 4, 1);
	printf("assignment:\n");
	fprint_uints(stdout, ic, n, 2, 1);
#endif
	
	/* Optimal transfer stage */
	for (unsigned int m = 0; m < max_iter; m++) {
		++(*niter);
#ifdef DEBUGGING_CODE
		printf("iter %d\n", m);
#endif
		for (i = 0; i < n; ++i) {
			indx++;

			/* do not move last member from a cluster */
			if (nclass[ic[i]] == 1)
				continue;

			/* In live cluster: try move to all clusters */
			if (i < aux_dp->live[ic[i]]) {
				aux_dp->index = aux_dp->all_index;
#ifdef DEBUGGING_CODE
				printf("%d is in live set %d\n", i, ic[i]);
#endif
				iter(d, auxd, seeds, K, n, p, nclass, ic, i, &indx);
			}
			/* Not in live cluster: try move only to live sets */
			else {
				aux_dp->index = aux_dp->some_index;
				len_index = 0;
				/* Find the live set */
				for (k = 0; k < K; ++k) {
					if (i < aux_dp->live[k])
						aux_dp->index[len_index++] = k;
					else if (k == ic[i])
						aux_dp->index[len_index++] = ic[i];
				}
#ifdef DEBUGGING_CODE
				printf("%d is in dead set %d\n", i, ic[i]);
				fprint_ints(stdout, aux_dp->live, len_index, 1, 1);
				
#endif
				iter(d, auxd, seeds, len_index, n, p, nclass, ic, i, &indx);
			}
			/* no transfers: algorithm has converged */
			if (indx == n)
				break;
		}
		for (k = 0; k < K; ++k)
			aux_dp->live[k] -= n;
		if (indx == n) {	/* [KSD] can move above preceding loop? */
			for (k = 0; k < K; ++k)
				for (j = 0; j < p; ++j)
					seeds[k][j] = aux_dp->last_cent[k][j];
			break;
		}
	}
	for (i = 0; i < n; ++i)
		for (j = 0; j < p; ++j)
			cost[ic[i]] +=
				dp->dmat[i][j] == dp->seeds[ic[i]][j]
				? aux_dp->prob_t[i][j]
				: aux_dp->prob_f[i][j];
/*
fprintf(stderr, "Update i:");
	for (i = 0; i < n; ++i)
fprintf(stderr, " %u", update_i[i]);
fprintf(stderr, "\nUpdate j:");
	for (j = 0; j < p; ++j)
fprintf(stderr, " %u", update_j[j]);
fprintf(stderr, "\n");
*/

	for (k = 0; k < K; ++k)
		sum += cost[k];
	return sum;
}/* cluster_hw_fast */


/**
 * Initialization for HW (fast).
 *
 * @param d		void pointer data structure: defined by user
 * @param seeds		seeds
 * @param nclass	count in clusters (1 x k)
 * @param ic		cluster assignments (1 x n)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 */
void
hw_fast_init_default (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic)
{
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int k, j, l;
	double max_llj;
	int max_id;

	/* initialize auxiliary variables: ekjn, vik, nclass, ic */
	init_auxiliary_variables_hw(d, auxd, seeds, K, n, p, nclass, ic);

	/* update haplotypes */
	for (k = 0; k < K; ++k) {
		for (j = 0; j < p; ++j) {
			max_llj = -INFINITY;
			max_id = -1;
			/* find largest e_kjN for each site */
			for (l = 0; l < dp->tot_n_categories; ++l)
				if (aux_dp->e_kjN[k][j][l] > max_llj) {
					max_llj = aux_dp->e_kjN[k][j][l];
					max_id = l;
				}
			seeds[k][j] = dp->categories[max_id];
			aux_dp->last_cent[k][j] = seeds[k][j];
		}
		aux_dp->live[k] = n + 1;	/* all clusters start to live */
	}

/*
	SETZERO_1ARRAY(nclass, K);
	for (unsigned int i = 0; i < n; ++i) {
		dp->last_assign[i] = ic[i];
		dp->idx_position[i] = nclass[ic[i]];
		dp->idx_in_cluster[ic[i]][nclass[ic[i]]++] = i;
	}
*/

} /* hw_fast_init_default */


/**
 * Optimal transfer stage.
 *
 * @param d		data pointer
 * @param seeds		centers
 * @param K		no. clusters to compare
 * @param n		no. observations
 * @param p		no. of coordinates
 * @param nclass	no. of observations in each cluster
 * @param ic		observation assignment
 * @param i		ith observation
 * @param indx		first 0, number of consecutive observations not transferred
 */
void
hw_fast_iter_default (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic,
	unsigned int i,
	unsigned int *indx)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = DEBUG_I;
#endif
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int k, k1, j, l;
	unsigned int last_assign = aux_dp->last_assign[i];
	double max_llj;
/* for debugging: remove when done
	double v_ik_diff;
	unsigned int i1;
*/
	int max_id;
	unsigned int cent_j_len;

/* for testing
if (!update_j)
	update_j = calloc(p, sizeof *update_j);

if (!update_i)
	update_i = calloc(n, sizeof *update_i);
*/
	
	double max_increase = INFINITY;
	/* Update the haplotypes and cost */
	for (k1 = 0; k1 < K; ++k1) {
		k = aux_dp->index[k1];
		aux_dp->llk_diff[k] = 0;
		cent_j_len = 0;
		for (j = 0; j < p; ++j) {
			/* [KSD] memory accesses of e_kjN and last_ekjn between iterations require swap */
			max_llj = -INFINITY;
			max_id = -1;
			
			for (l = 0; l < dp->tot_n_categories; ++l) {
				/* [KSD] Memory accesses in here are efficient */
				aux_dp->last_ekjn[k][j][l] = aux_dp->e_kjN[k][j][l];

				/* move to k */
				if (k != last_assign) {
					if (dp->dmat[i][j] != dp->categories[l])
						aux_dp->last_ekjn[k][j][l] += aux_dp->prob_f[i][j];
					else
						aux_dp->last_ekjn[k][j][l] += aux_dp->prob_t[i][j];
				}
				/* move from k */
				else {
					if (dp->dmat[i][j] != dp->categories[l])
						aux_dp->last_ekjn[k][j][l] -= aux_dp->prob_f[i][j];
					else
						aux_dp->last_ekjn[k][j][l] -= aux_dp->prob_t[i][j];
				}
				/* find the largest e_kjN[k][j][l] for each site*/
				if (aux_dp->last_ekjn[k][j][l] > max_llj) {
					max_llj = aux_dp->last_ekjn[k][j][l];
					max_id = l;
				}
			}
			seeds[k][j] = dp->categories[max_id]; //Compute new Hk
			
			/* adjust llk_diff for current observation */
			if (k != last_assign) {
				aux_dp->llk_diff[k] -= dp->dmat[i][j] == seeds[k][j]
					? aux_dp->prob_t[i][j]
					: aux_dp->prob_f[i][j];
				
			}
			else {
				aux_dp->llk_diff[k] -=
					dp->dmat[i][j] == aux_dp->last_cent[k][j]
						? aux_dp->prob_t[i][j]
						: aux_dp->prob_f[i][j];
			}

			if (seeds[k][j] != aux_dp->last_cent[k][j]) {
//update_j[j]++;
				/* Only consider the i1 in the current cluster k */
				if (k != last_assign)
					aux_dp->llk_diff[k] -=
						aux_dp->e_kjN[k][j][seeds[k][j]]
						- aux_dp->e_kjN[k][j][aux_dp->last_cent[k][j]];
				else
					aux_dp->llk_diff[k] +=
						aux_dp->last_ekjn[k][j][seeds[k][j]]
						- aux_dp->last_ekjn[k][j][aux_dp->last_cent[k][j]];
				aux_dp->changed_center[k][cent_j_len++] = j;
/*
double tmp = 0;
				for (i1 = 0; i1 < dp->cluster_size[dp->index[k]]; ++i1) {
					unsigned int idx = dp->idx_in_cluster[dp->index[k]][i1];
					if (idx != i) {
						v_ik_diff = dp->prob_t[idx][j] - dp->prob_f[idx][j];
						if (dp->dmat[idx][j] == seeds[dp->index[k]][j])
							//dp->llk_diff[dp->index[k]] += v_ik_diff;
							tmp += v_ik_diff;
						else if (dp->dmat[idx][j] == dp->last_cent[dp->index[k]][j])
							//dp->llk_diff[dp->index[k]] -= v_ik_diff;
							tmp -= v_ik_diff;
					}
				}
debug_msg(DEBUG_I, DEBUG_I, "%u (%u) %u (%u): %f %f %f\n", dp->index[k], dp->cluster_size[dp->index[k]], j, i, dp->index[k] != dp->last_assign[i] ? dp->e_kjN[dp->index[k]][j][seeds[dp->index[k]][j]] - dp->e_kjN[dp->index[k]][j][dp->last_cent[dp->index[k]][j]] : dp->last_ekjn[dp->index[k]][j][seeds[dp->index[k]][j]] - dp->last_ekjn[dp->index[k]][j][dp->last_cent[dp->index[k]][j]], tmp, (dp->index[k] != dp->last_assign[i] ? dp->e_kjN[dp->index[k]][j][seeds[dp->index[k]][j]] - dp->e_kjN[dp->index[k]][j][dp->last_cent[dp->index[k]][j]] : dp->last_ekjn[dp->index[k]][j][seeds[dp->index[k]][j]] - dp->last_ekjn[dp->index[k]][j][dp->last_cent[dp->index[k]][j]]) - tmp);
if (fabs((dp->index[k] != dp->last_assign[i] ? dp->e_kjN[dp->index[k]][j][seeds[dp->index[k]][j]] - dp->e_kjN[dp->index[k]][j][dp->last_cent[dp->index[k]][j]] : dp->last_ekjn[dp->index[k]][j][seeds[dp->index[k]][j]] - dp->last_ekjn[dp->index[k]][j][dp->last_cent[dp->index[k]][j]]) - tmp) > 1e-6)
exit(0);
*/
			}
		}
		/* only need to compare the changed log likelihodod */
/*
		if (k != last_assign)
			dp->llk_diff[k] = -dp->llk_diff[k];
*/
		if (aux_dp->llk_diff[k] < max_increase) {
			max_increase = aux_dp->llk_diff[k];
			ic[i] = k;
		}
		aux_dp->cent_j_len[k] = cent_j_len;
	}

	/* cluster assignment has changed */
	if (last_assign != ic[i]) {
//update_i[i]++;
		*indx = 0;
		move_observation(d, auxd, seeds, nclass, n, i, last_assign, ic[i]);
	}
} /* hw_fast_iter_default */

static void
move_observation (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int *nclass,
	unsigned int n,
	unsigned int i,
	unsigned int from,
	unsigned int to
	)
{
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int j, l;
	unsigned n_coord = dp->n_coordinates;
	unsigned n_categories = dp->tot_n_categories;
	unsigned int n_changed_j_to = aux_dp->cent_j_len[to];
	unsigned int n_changed_j_from = aux_dp->cent_j_len[from];

	aux_dp->idx_in_cluster[to][nclass[to]] = i;
	nclass[to]++;
	nclass[from]--;
		
	aux_dp->idx_in_cluster[from][aux_dp->idx_position[i]]
				= aux_dp->idx_in_cluster[from][nclass[from]];
	aux_dp->idx_position[aux_dp->idx_in_cluster[from][nclass[from]]]
							= aux_dp->idx_position[i];
	aux_dp->idx_position[i] = nclass[to] - 1;

	aux_dp->live[to] = aux_dp->live[from] = n + i + 1;

	/* Store the current center for the following comparision */
	for (j = 0; j < n_changed_j_to; ++j)
		aux_dp->last_cent[to][aux_dp->changed_center[to][j]]
			= seeds[to][aux_dp->changed_center[to][j]];
	for (j = 0; j < n_changed_j_from; ++j)
		aux_dp->last_cent[from][aux_dp->changed_center[from][j]]
			= seeds[from][aux_dp->changed_center[from][j]];
	for (j = 0; j < n_coord; ++j)
		for (l = 0; l < n_categories; ++l) {
			aux_dp->e_kjN[to][j][l] = aux_dp->last_ekjn[to][j][l];
			aux_dp->e_kjN[from][j][l] = aux_dp->last_ekjn[from][j][l];
		}
	/* Store the last assignment */
	aux_dp->last_assign[i] = to;
} /* move_observation */


/**
 * Initialize log read probabilities for each cluster v_{ik} and auxiliary
 * variables ekjn, as well as cluster assignments and cluster sizes.
 *
 * @param d		void pointer to data structure: defined by user
 * @param seeds		seeds
 * @param K		number of clusters
 * @param n		number of observations
 * @param p		number of coordinates
 * @param nclass	count in clusters (1 x k)
 * @param ic1		cluster assignments (1 x n)
 */
static void
init_auxiliary_variables (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic1)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;
#endif
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int i, k, j, l;

	/* initialize */
	SETZERO_1ARRAY(nclass, K);

	for (k = 0; k < K; ++k)
		for (j = 0; j < p; ++j)
			for (l = 0; l < dp->tot_n_categories; ++l)
				aux_dp->e_kjN[k][j][l] = 0;
	
	for (i = 0; i < n; ++i) {
		/* identify closest haplotype */
		double maxll = -INFINITY;
		int max_id = -1;
		for (k = 0; k < K; ++k) {    /* compare to each haplotype */
			aux_dp->v_ik[i][k] = 0;
			for (j = 0; j < p; ++j)
				aux_dp->v_ik[i][k] +=
					dp->dmat[i][j] == seeds[k][j]
					? aux_dp->prob_t[i][j]
					: aux_dp->prob_f[i][j];
#ifdef DEBUGGING_CODE
			debug_msg(DEBUG_II, fxn_debug, "v[%u][%u]: %f\n",
						i, k, aux_dp->v_ik[i][k]);
#endif
			/* find the biggest v_ik for each read */
			if (aux_dp->v_ik[i][k] > maxll) {
				maxll = aux_dp->v_ik[i][k];
				max_id = k;
			}
		}
		nclass[max_id]++;
		ic1[i] = max_id;
		
		/* update log probability of each cluster, pos., hap. nuc. */
		for (j = 0; j < p; ++j)
			for (l = 0; l < dp->tot_n_categories; ++l)
				aux_dp->e_kjN[ic1[i]][j][l] +=
					dp->dmat[i][j] == dp->categories[l]
					? aux_dp->prob_t[i][j] : aux_dp->prob_f[i][j];
	}
} /* init_auxiliary_variables */

/**
 * Initialize auxiliary variables used by efficient hw algorithm, as well as
 * cluster assignments and cluster sizes.
 *
 * @param d		void pointer to data structure: defined by user
 * @param seeds		seeds
 * @param K		number of clusters
 * @param n		number of observations
 * @param p		number of coordinates
 * @param nclass	count in clusters (1 x k)
 * @param ic1		cluster assignments (1 x n)
 */
static void
init_auxiliary_variables_hw (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p,
	unsigned int *nclass,
	unsigned int *ic1)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;
#endif
	data *dp = (data *) d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int i, k, j, l;

	SETZERO_1ARRAY(nclass, K);

	/* initialize e_kjN */
	for (k = 0; k < K; ++k)
		for (j = 0; j < p; ++j)
			for (l = 0; l < dp->tot_n_categories; ++l)
				aux_dp->e_kjN[k][j][l] = 0;
	
	for (i = 0; i < n; ++i) {
		/* identify closest haplotype */
		double maxll = -INFINITY;
		int max_id = -1;
		for (k = 0; k < K; ++k) {    /* compare to each haplotype */
			aux_dp->v_ik[i][k] = 0;
			for (j = 0; j < p; ++j)
				aux_dp->v_ik[i][k] +=
					dp->dmat[i][j] == seeds[k][j]
					? aux_dp->prob_t[i][j]
					: aux_dp->prob_f[i][j];
#ifdef DEBUGGING_CODE
			debug_msg(DEBUG_II, fxn_debug, "v[%u][%u]: %f\n",
						i, k, aux_dp->v_ik[i][k]);
#endif
			/* find the biggest v_ik for each read */
			if (aux_dp->v_ik[i][k] > maxll) {
				maxll = aux_dp->v_ik[i][k];
				max_id = k;
			}
		}

		/* set up assignment and various indices */
		ic1[i] = max_id;
		aux_dp->last_assign[i] = max_id;
		aux_dp->idx_position[i] = nclass[max_id];
		aux_dp->idx_in_cluster[max_id][nclass[max_id]] = i;
		nclass[max_id]++;
		
		/* update log probability of each cluster, pos., hap. nuc. */
		for (j = 0; j < p; ++j)
			for (l = 0; l < dp->tot_n_categories; ++l)
				aux_dp->e_kjN[ic1[i]][j][l] +=
					dp->dmat[i][j] == dp->categories[l]
					? aux_dp->prob_t[i][j] : aux_dp->prob_f[i][j];
	}

} /* init_auxiliary_variables_hw */


/**
 * Initialize the haplotypes and vik.
 *
 * @param d	void pointer data structure: defined by user
 * @param seeds	seeds
 * @param K	number of clusters
 * @param n	number of observations
 * @param p	number of coordinates
 */
static void
init_haplotypes (
	void *d,
	void *auxd,
	data_t **seeds,
	unsigned int K, unsigned int n, unsigned int p)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
#endif
	data *dp = (data *)d;
	auxvar *aux_dp = (auxvar *) auxd;
	unsigned int i, k, j, l;
	double v_ik_diff;

	/* compute most likely haplotype nucleotides and update vik */
	for (k = 0; k < K; ++k)
		for (j = 0; j < p; ++j) {
			aux_dp->last_cent[k][j] = seeds[k][j];
			double max_llj = -INFINITY;
			int max_id = -1;
			for (l = 0; l < dp->tot_n_categories; ++l)
				/* find largest e_kjN for each site */
				if (aux_dp->e_kjN[k][j][l] > max_llj) {
					max_llj = aux_dp->e_kjN[k][j][l];
					max_id = l;
				}
			seeds[k][j] = dp->categories[max_id];

			/* update vik */
			if (seeds[k][j] != aux_dp->last_cent[k][j]) {
				for (i = 0; i < n; ++i) {
					v_ik_diff = aux_dp->prob_t[i][j]
						- aux_dp->prob_f[i][j];
					if (dp->dmat[i][j] == seeds[k][j])
						aux_dp->v_ik[i][k] += v_ik_diff;
					else if (dp->dmat[i][j]
						== aux_dp->last_cent[k][j])
						aux_dp->v_ik[i][k] -= v_ik_diff;
				}
				aux_dp->last_cent[k][j] = seeds[k][j];
			}
		}
} /* init_haplotypes */

/**
 * Find the unique sequence with the minimum llk under current assignment, make it as the next seed.
 *
 * @param d data	structure
 * @param length	no. of sequence chosen in that hash table
 * @param abundance	abundance of each seq
 * @param unique_seq	index of each seq
 * @param seq_arr	abundance idx of each seq
 */
void
seed_llk (
	void *d,
	void *auxd,
	unsigned int length,
	unsigned int *abundance,
	size_t *unique_seq,
	size_t ** seq_arr)
{
	unsigned int i, j, k, l;
	unsigned int keep_going;
	data *dp = (data *)d;
	auxvar *aux_dp = (auxvar *) auxd;

	SETZERO_1ARRAY(dp->llk_seq, length);
	
	for (i = 0; i < length; ++i)
		for (j = 0; j < abundance[i]; ++j)
			if (dp->cluster_id[unique_seq[i]]
				== dp->cluster_id[seq_arr[i][j]])
				dp->llk_seq[i] += aux_dp->v_ik[seq_arr[i][j]][dp->cluster_id[unique_seq[i]]];
	
	double min = 0;
	size_t min_id = -1;
	
	/* Find the minimum llk, k starts from 1 to exclude the most one that has been choosed before */
	for (k = 1; k < length; ++k) {
		keep_going = 1;
		if (dp->seed_count != 0)
			for (l = 0; l < dp->seed_count; ++l)
				if (unique_seq[k] == dp->min_index[l])
					keep_going = 0;
		if (keep_going == 1)
			if (dp->llk_seq[k] < min) {
				min = dp->llk_seq[k];
				min_id = unique_seq[k];
			}
	}
	dp->min_index[dp->seed_count] = min_id;
	dp->seed_count++;

	COPY_1ARRAY(dp->seeds[dp->seed_count], dp->dmat[min_id],
						dp->n_coordinates);
}/* seed_llk */


/**
 * Human-friendly string defining initialization method.
 *
 * @param init_method	id of initialization method
 * @return		constant character string
 */
char const *khaplotype_init_method(int init_method)
{
	if (init_method < KMODES_INIT_NUMBER_METHODS)
		return kcluster_init_method(init_method);
	else if (init_method == KHAPLOTYPE_INIT_FILTER_DATA)
		return "initialize by filtering the data with hash and expected errors";
	return "unknown initialization method";
} /* khaplotype_init_method */


/**
 * Returns human-friendly description of a k-haplotype error.
 *
 * @param err integer representation of error; something returned by
 *	    k-haplotype methods or initialiation methods.
 * @return descriptive string explaining error or NULL if invalid error
 */
const char *khaplotype_error(int err) {
	if (err == KHAPLOTYPE_EXCEED_ITER_WARNING)
		return "Exceeded maximum iterations";
	else if (err == KHAPLOTYPE_NULL_CLUSTER_ERROR)
		return "Inferred null cluster";
	else if (err == KHAPLOTYPE_CALLER_INPUT_ERROR)
		return "Requested K=0 or K>n, more clusters than observations";
	else if (err == KHAPLOTYPE_MEMORY_ERROR)
		return "Memory allocation error";
	else if (err == KHAPLOTYPE_INVALID_INITIALIZATION_METHOD)
		return "Unknown initialization method requested";
	else if (err == KHAPLOTYPE_INTERNAL_ERROR)
		return "Internal error";
	else
		return NULL;
} /* khaplotype_error */

