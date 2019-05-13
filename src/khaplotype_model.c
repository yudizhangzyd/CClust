#include "khaplotype_model.h"
#include "math.h"

/**
 * Create model object.
 *
 * @param mod	model object to create
 * @param dat	pointer to data object
 * @return	error status
 */
int make_model(model **mod, data *dat)
{
	int err = NO_ERROR;
	model *rm;
	*mod = malloc(sizeof **mod);
	
	if (*mod == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model");
	rm = *mod;
	
	rm->haplotypes = NULL;
	rm->haplotypes = malloc(dat->max_read_length * (dat->seed_count + 1)
				* sizeof *rm->haplotypes);
	
	if (!rm->haplotypes)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"model::haplotypes");
	
	rm->samp_prob = NULL;
	rm->samp_prob = malloc(dat->n_observations * sizeof(*rm->samp_prob));
	
	rm->JC_ll = -INFINITY;
	rm->best_ll = -INFINITY;
	
	rm->aic = INFINITY;
	rm->bic = INFINITY;
	rm->previous_bic = INFINITY;
	rm->percent_change = 0.0;
	
	/* include the 1st seed */
	rm->n_param = (dat->seed_count + 1) * dat->max_read_length;	/* haplotypes */
	//	+ dat->seed_count;				/* pi, [QUESTION: How to choose pars? Should I include the first haplotype I choose?] */
	
	rm->distance = NULL;
	rm->est_ancestor = NULL;
	
	rm->distance = malloc((dat->seed_count + 1) * sizeof *rm->distance);
	rm->est_ancestor = malloc(dat->max_read_length
				  * sizeof * rm->est_ancestor);
	
	if (!rm->distance || !rm->est_ancestor)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"model::JC69");
	return err;
}/* make_model */

/**
 * Reallocate model struct for a different K or different sample with different size.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @return	error status
 */
int realloc_model(model *mod, data *dat)
{
	unsigned char *haplotypes = realloc(mod->haplotypes,
					    dat->max_read_length * (dat->seed_count + 1) * sizeof *mod->haplotypes);
	
	if (!haplotypes)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"realloc.model.haplotypes");
	
	for (unsigned int i = 0; i < dat->seed_count + 1; ++i)
		for (unsigned int j = 0; j < dat->max_read_length; ++j)
			haplotypes[i*dat->max_read_length + j] = dat->seeds[i][j];

	mod->haplotypes = haplotypes;
	
	mod->n_param = (dat->seed_count + 1) * dat->max_read_length;	/* haplotypes */
	//	+ dat->seed_count;				/* pi, [QUESTION: How to choose pars? Should I include the first haplotype I choose?] */
	mod->previous_bic = mod->bic;
	mod->aic = INFINITY;
	mod->bic = INFINITY;
	
	mod->distance = realloc(mod->distance, (dat->seed_count + 1) * sizeof *mod->distance);
	
	mmessage(INFO_MSG, NO_ERROR, "Number of parameters: %u\n", mod->n_param);
	
	return NO_ERROR;
}/* realloc_model */

/**
 * Delete model object.
 *
 * @param mod	pointer to model object to delete
 */
void free_model(model *mod)
{
	if (!mod)
		return;
	if (mod->haplotypes) free(mod->haplotypes);
	if (mod->samp_prob) free(mod->samp_prob);
	if (mod->distance) free(mod->distance);
	if (mod->est_ancestor) free(mod->est_ancestor);
	
	free(mod);
	mod = NULL;
} /* free_model */

/**
 * Calculate aic and bic modified by approximate JC69 hierarchical model on
 * haplotypes.
 *
 * @param hap			haplotypes
 * @param est_anc		ancester sequence
 * @param distance		distance from haplotypes to ancestor sequence
 * @param best_ll		current log likelihood from data
 * @param K			number of haplotypes
 * @oaram JC_ll			log likelihood from JC69 model
 * @param n_aic			pointer to aic, value updated
 * @param n_bic			pointer to bic, value updated
 * @param n_param		number of parameters in current model
 * @param max_read_length	length of haplotypes
 * @param sample_size		sample size
 *
 * return			err status
 **/
int modified_ic(unsigned char *hap, unsigned char *est_anc, double *distance,
		double best_ll, unsigned int K, double *JC_ll, double *n_aic,
		double *n_bic, unsigned int n_param, unsigned int max_read_length,
		size_t sample_size)
{
	int param_change = 0;
	
	m_JC69(hap, est_anc, distance, K, max_read_length);
	*JC_ll = e_JC69(hap, est_anc, distance, K, max_read_length);
	
	/* K branch lengths, ancestral haplotype, but no haplotypes estimated */
	param_change = K - max_read_length * (K - 1);
	
	*n_aic = aic(best_ll + *JC_ll, n_param + param_change);
	*n_bic = bic(best_ll, n_param, sample_size);
	
	return NO_ERROR;
}/* modified_ic */

/**
 * get MMEs of all parameters in the JC69 model.
 *
 * @param hap	pointer to haplotypes
 * @param anc	ancestor haplotype, to be calculated (tbc)
 * @param dist	expected no. changes/site b/w ancestor & haplotype, tbc
 * @param K	number of haplotypes
 * @param len	length of reads
 * @return err
 **/
int m_JC69(unsigned char * hap, unsigned char * anc, double *dist,
	   unsigned int K, unsigned int len)
{
	
	int err = NO_ERROR;
	unsigned int count[NUM_NUCLEOTIDES];
	unsigned int max_count;
	
	/* most common nucleotide across haplotypes is the estimated ancestor */
	for (unsigned int j = 0; j < len; j ++){
		for (unsigned char n = 0; n < NUM_NUCLEOTIDES; n++)
			count[n] = 0;
		for (unsigned int k = 0; k < K; k++)
			count[hap[k * len + j]]++;
		max_count = 0;
		for (unsigned char n = 0; n < NUM_NUCLEOTIDES; ++n)
			if (count[n] > max_count) {
				max_count = count[n];
				anc[j] = n;
			}
	}
	
	/* [KSD, BUG] Was the expected number of OBSERVED changes per site. */
	/* estimate the expected number of changes per site of all haplotypes */
	for (unsigned int k = 0; k < K; k++) {
		if (1) {	/* bug-free version */
			double tmp = (double) hamming_char_dis( (char *) &hap[k*len],
							       (char *) anc, (size_t) len) / len;
			dist[k] = -0.75 * log(1 - tmp / 0.75);
		} else {	/* buggy version */
			dist[k] = (double) hamming_char_dis( (char *) &hap[k*len], (char *) anc, (size_t) len )/len;
		}
	}
	
	return err;
}/* m_JC69 */

/**
 * Calculate the log likelihood of all K haplotypes under JC69 model.
 *
 * @param hap   haplotype sequences
 * @param anc   ancestor sequence
 * @param dist	expected no. changes/site for each haplotype
 * @param K	number of haplotypes
 * @param len	length of the haplotypes
 * @return	err status
 **/
double e_JC69(unsigned char * hap, unsigned char * anc, double *dist,
	      unsigned int K, unsigned int len)
{
	double ll = 0;
	
	for (unsigned int k = 0; k < K; k++)
		for (unsigned int j = 0; j < len; j++)
			if (anc[j] == hap[k * len + j])
				ll += log(0.25 + 0.75 * exp(-dist[k] / 0.75));
			else
				ll += log(0.25 - 0.25 * exp(-dist[k] / 0.75));
	
	return ll;
}/* e_JC69 */

int sample_from_categorical_distribution(data *dat, model *mod)
{
	double r, cdf;
	r = rand() / (RAND_MAX + 1.);
	
	mod->ind = 0;
	cdf = mod->samp_prob[mod->ind];
	while (r > cdf) {
		cdf += mod->samp_prob[++mod->ind];
		if (mod->ind == dat->n_observations - 1)
			break;
	}
	return mod->ind;
}

/**
 * check how many true modes the algorithms find

 @param optimal_seq modes get by algorithm
 @param true_seq true modes
 @param K no. of modes
 @param p coordinates of modes
 @return no. of match
 */
int detect_true (data_t **optimal_seq, data_t **true_seq, unsigned int K, unsigned int p) {
	unsigned int i, k, k1, j, length = 0;
	unsigned int *match = NULL;
	MAKE_1ARRAY(match, K);
	for (k = 0; k < K; ++k)
		for (k1 = 0; k1 < K; ++k1){
			int equal = 1;
			for (j = 0; j < p; ++j)
				if (optimal_seq[k][j] != true_seq[k1][j])
					equal = 0;
			if (equal)
				match[length++] = k1;
		}
	
	for(i = 0; i < length; ++i)
		for(j = i + 1; j < length;)
		{
			if(match[i] == match[j])
			{
				for(k = j; k < length - 1; ++k)
					match[k] = match[k + 1];
				--length;
			}
			else
				++j;
		}
	free(match);
	return length;
}
