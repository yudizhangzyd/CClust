/**
 * @file khaplotype_initialize.c
 *
 * Initialization code for k-haplotype method.  The methods defined here are
 * aware the data are fastq data.
 */

#include "khaplotype.h"
#include "khaplotype_initialize.h"

static int khaplotype_init_filter_data(data *dat, auxvar *aux_dat, options *opt, model *mod);

/**
 * Initialize seeds.  These methods work for k-modes, but data-specific methods
 * can be provided in the call-back function.
 *
 * @param dat	data structure
 * @param opt	option structure
 * @return	err status
 */
int initialize(data *dat, auxvar *aux_dat, options *opt, model *mod)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
#endif
	int err = NO_ERROR;
	data_t **seeds = dat->seeds;
	unsigned int *seed_idx = dat->seed_idx;

	if (opt->n_inner_init > 1)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
			"Not set up to do nested initialization.\n");

	for (unsigned int j = 0; j < opt->n_inner_init; ++j) {
		if (opt->init_method
		    == KHAPLOTYPE_INIT_FILTER_DATA)
			khaplotype_init_filter_data(dat, aux_dat, opt, mod);
		else if (opt->init_method < KMODES_INIT_NUMBER_METHODS)
			kcluster_init(dat->dmat, dat->n_observations,
				      dat->n_coordinates, opt->K,
				      opt->n_seed_set, seeds, seed_idx,
				      opt->init_method);
		else
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Invalid initialization method");
	}

	return err;
} /* initialize */


/**
 * Use hash and expected mean errors to filter the data and do the initialization.
 *
 *
 * @param dat	data structure
 * @param opt	option structure
 * @param mod	model structure
 * @return	error status
 */
static int
khaplotype_init_filter_data (
	data *dat,
	auxvar *aux_dat,
	options *opt,
	model *mod)
{
	int err = NO_ERROR;
	int first = 1;
	dat->seed_count = 0;
	dat->use_ini = 1;
	double sum_prob = 0;

	/* use the most abundant haplotype as the first seed */
	COPY_1ARRAY(dat->seeds[0], dat->dmat[dat->seq_count->idx], dat->n_coordinates);
	SETZERO_1ARRAY(dat->cluster_size, dat->seed_count + 1);
	SETZERO_2ARRAY(aux_dat->v_ik, opt->K);
	fastq_lloyds_efficient_step1(dat, aux_dat, dat->seeds, dat->seed_count + 1, dat->n_observations, dat->n_coordinates, dat->cluster_size, dat->cluster_id, 1);
	
	/* compute the BIC for the 1st seed */
	mod->best_ll = compute_cost(dat, dat->cluster_id, dat->criterion, dat->seed_count + 1, dat->n_observations);
	mod->bic = bic(mod->best_ll, mod->n_param,
		       dat->n_observations);
	
	/* add the haplotypes by identifying the least llk */
	if (opt->filter_haplotypes)
		MAKE_1ARRAY(dat->llk_seq, dat->n_selected_sequences);
	else
		MAKE_1ARRAY(dat->llk_seq, dat->hash_leng_exk);
	
	for (unsigned int k1 = 0; k1 < opt->K - 1; ++k1) {
		
		if (opt->filter_haplotypes)
			seed_llk(dat, aux_dat, dat->n_selected_sequences, dat->sel_uniq_seq_count, dat->sel_uniq_seq_idx, dat->sel_idx_arr);
		else
			seed_llk(dat, aux_dat, dat->hash_leng_exk, dat->seq_abundance, dat->hash_idx_exk, dat->hash_idx_exk_arr);
		
		SETZERO_1ARRAY(dat->cluster_size, dat->seed_count +1);
		SETZERO_2ARRAY(aux_dat->v_ik, opt->K);
		fastq_lloyds_efficient_step1(dat, aux_dat, dat->seeds, dat->seed_count +1, dat->n_observations, dat->n_coordinates, dat->cluster_size, dat->cluster_id, 1);
		
		/* Compute BIC to decide a cut-off value for randomnization */
		if ((err = realloc_model(mod, dat)))
			return err;
		mod->best_ll = compute_cost(dat, dat->cluster_id, dat->criterion, dat->seed_count + 1, dat->n_observations);
		modified_ic(mod->haplotypes, mod->est_ancestor,
			    mod->distance, mod->best_ll, dat->seed_count + 1,
			    &mod->JC_ll, &mod->aic, &mod->bic,
			    mod->n_param, dat->max_read_length,
			    dat->n_observations);
		mod->percent_change = (mod->previous_bic - mod->bic) / mod->previous_bic;
		
		if(first)
			if (mod->percent_change < opt->cut_k) {
				opt->true_K = dat->seed_count + 1;
				first = 0;
			}
		
		if (dat->seed_count == opt->K - 1)
			break;
	}
	
	/* Now decide the cut-off value depending on BIC, despite two models are indistinguishable by BIC criterion if the difference < 2, now we only have two situations. 1st, cluster K > true K, randomnization plays a role, otherwise choose all of them  */
	if (opt->K <= opt->true_K) {
				for (unsigned int i = 0; i < opt->K; ++i)
					COPY_1ARRAY(dat->ini_seeds[i], dat->seeds[i], dat->n_coordinates);
			}
	else{
		for (unsigned int i = 0; i < opt->true_K; ++i)
			COPY_1ARRAY(dat->ini_seeds[i], dat->seeds[i], dat->n_coordinates);
		
		unsigned int start = 0 ,length_id = 0, id = 0;
		unsigned int *previous_id = NULL;
		previous_id = malloc((length_id + 1) * sizeof(*previous_id));
		for (unsigned int k1 = 0; k1 < opt->K - opt->true_K; ++k1) {
			
			/* store the initial sampling probablity */
			for (unsigned int l = 0; l < dat->n_observations; ++l) {
				mod->samp_prob[l] = aux_dat->v_ik[l][dat->cluster_id[l]];
				sum_prob += mod->samp_prob[l];
			}
			
			for (unsigned int l = 0; l < dat->n_observations; ++l)
				mod->samp_prob[l] = mod->samp_prob[l]/sum_prob;
			
			for (unsigned int m = 0; m < dat->max_iter; ++m) {
				id = sample_from_categorical_distribution(dat, mod);
				unsigned int keep_going = 1;
				
				if (start)
					for (unsigned int i = 0; i < length_id; ++i)
						if (id == previous_id[i])
							keep_going = 0;
				if (id == dat->hash_idx_exk[0])
					keep_going = 0;
				else
					for (unsigned int j = 0; j < opt->true_K; ++j)
						if (id == dat->min_index[j])
							keep_going = 0;
				if (keep_going) {
					start = 1;
					COPY_1ARRAY(dat->ini_seeds[k1 + opt->true_K], dat->dmat[id], dat->n_coordinates);
					break;
				}
			}
			previous_id[++length_id] = id;
			previous_id = realloc(previous_id, (length_id + 1) * sizeof(*previous_id));
			
			SETZERO_1ARRAY(dat->cluster_size, k1 + opt->true_K + 1);
			SETZERO_2ARRAY(aux_dat->v_ik, opt->K);
			fastq_lloyds_efficient_step1(dat, aux_dat, dat->seeds, k1 + opt->true_K + 1, dat->n_observations, dat->n_coordinates, dat->cluster_size, dat->cluster_id, 1);
		}
		free(previous_id);
	}
	return err;
}/* khaplotype_init_filter_data */
