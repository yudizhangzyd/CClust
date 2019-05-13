/**
 * @file haplotype.c
 * @author Yudi Zhang
 *
 * [TODO] Improve use of hash tables.  Try to reuse the same code?
 */

#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>

#include "khaplotype.h"
#include "khaplotype_data.h"
#include "khaplotype_initialize.h"
#include "cluster_lloyds.h"
#include "array.h"
#include "error.h"
#include "io.h"
#include "math.h"

/**
 * Sample n from N objects without replacemnt.
 *
 * @param N	number of objects to sample from
 * @param n	number of objects to sample
 * @param idx	0-based indices of n selected objects from set {0,1,...,N-1}
 */
void sample_better(unsigned int N, unsigned int n, unsigned int *idx)
{
	unsigned int t = 0, m = 0, lim, i;

	while (m < n) {
		/* note integer division */
		lim = (N - t) * (RAND_MAX / (N - t));
		do {
			i = rand();
		} while (i >= lim);
		i %= (N - t);
		if (i < n - m) {
//			unsigned int keep_going = 1;
//			if (t == dat->hash_idx_exk[0])
//				keep_going = 0;
//			else
//				for (unsigned int j = 0; j < dat->seed_count; ++j)
//					if (t == dat->min_index[j])
//						keep_going = 0;
//			if (keep_going)
				idx[m++] = t++;
		}
		else
			++t;
	}
} /* sample_better */


int sync_state_aux_var(data *dat, auxvar *aux, options *opt)
{
	/* initialize the auxiliary variables */
	if (opt->run_with_quals) {
		MAKE_2ARRAY(aux->prob_f, dat->n_observations, dat->max_read_length);
		MAKE_2ARRAY(aux->prob_t, dat->n_observations, dat->max_read_length);
		compute_pij(dat, aux->prob_t, aux->prob_f);
	}
	else {
		MAKE_2ARRAY(aux->prob_f, dat->n_observations, dat->n_coordinates);
		MAKE_2ARRAY(aux->prob_t, dat->n_observations, dat->n_coordinates);
		compute_pij(dat, aux->prob_t, aux->prob_f);
	}
	
	MAKE_2ARRAY(aux->v_ik, dat->n_observations, opt->K);
	MAKE_3ARRAY(aux->e_kjN, opt->K, dat->n_coordinates, dat->tot_n_categories);
	MAKE_1ARRAY(aux->last_assign, dat->n_observations);
	MAKE_2ARRAY(aux->last_cent, opt->K, dat->n_coordinates);
	MAKE_1ARRAY(aux->last_cost, opt->K);
	MAKE_2ARRAY(aux->last_vik, dat->n_observations, opt->K);
	MAKE_3ARRAY(aux->last_ekjn, opt->K, dat->n_coordinates, dat->tot_n_categories);
	MAKE_1ARRAY(aux->live, opt->K);
	MAKE_1ARRAY(aux->some_index, opt->K);
	MAKE_1ARRAY(aux->all_index, opt->K);
	MAKE_1ARRAY(aux->llk_diff, opt->K);
	MAKE_2ARRAY(aux->idx_in_cluster, opt->K, dat->n_observations);
	MAKE_1ARRAY(aux->index_in_cluster, dat->n_observations);
	MAKE_1ARRAY(aux->idx_position, dat->n_observations);
	MAKE_2ARRAY(aux->changed_center, opt->K, dat->n_coordinates);
	MAKE_1ARRAY(aux->cent_j_len, opt->K);
	
	aux->index = NULL;
	for (unsigned int k = 0; k < opt->K; ++k)
		aux->all_index[k] = k;
	
	return NO_ERROR;
} /* sync_state_aux_var */

int make_aux_var(auxvar **aux)
{
	*aux = malloc(sizeof **aux);
	
	if (*aux == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "aux variables object");
	(*aux)->prob_f = NULL;
	(*aux)->prob_t = NULL;
	(*aux)->v_ik = NULL;
	(*aux)->e_kjN = NULL;
	(*aux)->last_cent = NULL;
	(*aux)->last_assign = NULL;
	(*aux)->last_cost = NULL;
	(*aux)->last_vik = NULL;
	(*aux)->last_ekjn = NULL;
	(*aux)->live = NULL;
	(*aux)->index = NULL;
	(*aux)->index_in_cluster = NULL;
	(*aux)->llk_diff = NULL;
	(*aux)->idx_in_cluster = NULL;
	(*aux)->idx_position = NULL;
	(*aux)->changed_center = NULL;
	(*aux)->cent_j_len = NULL;
	return NO_ERROR;
}/* make_aux_var */

void free_aux_var(auxvar *aux)
{
	if (!aux) return;
	
	FREE_2ARRAY(aux->prob_f);
	FREE_2ARRAY(aux->prob_t);
	FREE_2ARRAY(aux->v_ik);
	FREE_3ARRAY(aux->e_kjN);
	FREE_1ARRAY(aux->last_assign);
	FREE_2ARRAY(aux->last_cent);
	FREE_1ARRAY(aux->last_cost);
	FREE_2ARRAY(aux->last_vik);
	FREE_3ARRAY(aux->last_ekjn);
	FREE_1ARRAY(aux->live);
	FREE_VECTOR(aux->some_index);
	FREE_VECTOR(aux->all_index);
	FREE_VECTOR(aux->index_in_cluster);
	FREE_VECTOR(aux->llk_diff);
	FREE_2ARRAY(aux->idx_in_cluster);
	FREE_1ARRAY(aux->idx_position);
	FREE_2ARRAY(aux->changed_center);
	FREE_1ARRAY(aux->cent_j_len);
}/* free_aux_var */

/**
 * Sync state of data object.  Allocate auxiliary variables.
 *
 * @param dat	pointer to data object
 * @param opt	pointer to option object
 * @return	error status
 */
int sync_state_khaplotype_data(data *dat, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;

	/* get coordinates from data */
	dat->n_observations = dat->fdata->n_reads;

	dat->max_read_length = dat->fdata->n_max_length;
	dat->min_read_length = dat->fdata->n_min_length;

	if (dat->max_read_length == dat->min_read_length)
		dat->n_coordinates = dat->min_read_length;

	dat->n_quality = dat->fdata->max_quality - dat->fdata->min_quality + 1;
	dat->lengths = malloc(dat->n_observations * sizeof *dat->lengths);

	if (!dat->lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data.lengths");

	if (dat->fdata->n_lengths) {
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "k-haplotypes "
				"does not work with variable length reads!");
		/*
		memcpy(dat->lengths, dat->fdata->n_lengths, dat->n_observations
			* sizeof *dat->lengths);
		*/
	} else {
		for (size_t i = 0; i < dat->n_observations; ++i)
			dat->lengths[i] = dat->max_read_length;
	}

	/* allocate the index array of reads */
	dat->read_idx = malloc(dat->fdata->n_reads * sizeof *dat->read_idx);

	if (!dat->read_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.read_idx");

	for (size_t i = 0; i < dat->fdata->n_reads; ++i)
		dat->read_idx[i] = i;

	/* make reads matrix */
	dat->dmat = malloc(dat->n_observations * sizeof *dat->dmat);

	if (!dat->dmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.dmat");

	unsigned char *rptr = dat->fdata->reads;
	for (size_t i = 0; i < dat->n_observations; ++i) {
		dat->dmat[i] = rptr;
		rptr += dat->lengths[i];
	}

	debug_msg(DEBUG_I, fxn_debug, "Allocated %dx(%d) sequence matrix\n",
				dat->n_observations, dat->max_read_length);

	/* make quals matrix */
	/* allocate short-cut pointers to quality sequences */
	dat->qmat = malloc(dat->n_observations * sizeof *dat->qmat);

	if (!dat->qmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.qmat");

	unsigned char *qptr = dat->fdata->quals;
	for (size_t i = 0; i < dat->n_observations; i++) {
		dat->qmat[i] = qptr;
		qptr += dat->lengths[i];
	}

	debug_msg(DEBUG_I, fxn_debug, "Allocated %dx(%d) quality matrix\n",
		dat->n_observations, dat->max_read_length);

	dat->categories = find_unique(dat->fdata->reads, dat->n_observations
				* dat->max_read_length, &dat->tot_n_categories);
	
	/* Allocate for writing results */
	dat->cluster_id = malloc(dat->n_observations * sizeof *dat->cluster_id);
	dat->best_cluster_id = malloc(dat->n_observations * sizeof *dat->best_cluster_id);

	if (dat->cluster_id == NULL || dat->best_cluster_id == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::cluster_id");

	dat->cluster_size = malloc(opt->K * sizeof *dat->cluster_size);
	dat->best_cluster_size = malloc(opt->K * sizeof *dat->best_cluster_size);

	if (dat->cluster_size == NULL || dat->best_cluster_size == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::cluster_size");

	dat->criterion = malloc(opt->K * sizeof *dat->criterion);
	dat->best_criterion = malloc(opt->K * sizeof *dat->best_criterion);

	if (dat->criterion == NULL || dat->best_criterion == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::criterion");

	if (opt->shuffle) {
		dat->obsn_idx = malloc(dat->n_observations
					* sizeof *dat->obsn_idx);
		dat->best_obsn_idx = malloc(dat->n_observations
					* sizeof *dat->obsn_idx);
		if (!dat->obsn_idx || !dat->best_obsn_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"data:obsn_idx");
		for (unsigned int i = 0; i < dat->n_observations; ++i)
			dat->obsn_idx[i] = i;
	}

	if (opt->sim_info_file) {
		opt->sim_cluster = malloc(dat->n_observations
					* sizeof *opt->sim_cluster);
		if (!opt->sim_cluster)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"options:sim_cluster");
		FILE *fp = fopen(opt->sim_info_file, "r");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
					opt->sim_info_file);
		if ((err = fread_uints(fp, opt->sim_cluster,
				dat->n_observations))) {
			free(opt->sim_cluster);
			return err;
		}
		fclose(fp);
		opt->sim_K = 0;
		for (unsigned int i = 0; i < dat->n_observations; ++i)
			if (opt->sim_cluster[i] > opt->sim_K)
				opt->sim_K = opt->sim_cluster[i];
		++opt->sim_K;
	}
	
//	if (opt->true_modes_file) {
//		MAKE_2ARRAY(opt->sim_modes, opt->K, dat->n_coordinates);
//		if ((err = read_fsa(opt->true_modes_file, opt->K, dat->n_coordinates, opt->sim_modes))) {
//			free(opt->sim_modes);
//			return err;
//		}
//	}
	
	return err;
}/* sync_state_khaplotype_data */


/**
 * Set up seeds

 * @param dat data
 * @param opt option
 * @return err status
 */
int make_seeds(data *dat, options *opt)
{

	if (!opt->K)
		return NO_ERROR;

	if (!opt->seed_idx) {
		dat->seed_idx = malloc(opt->K * sizeof *dat->seed_idx);
		if (!dat->seed_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"data::seeds_idx");
	} else {
		dat->seed_idx = opt->seed_idx;
	}

/* [BUG, TODO] data::ini_seeds is not being used correctly here, esp.
 * when opt->n_inner_init > 1, in which case this is a memory leak.
 */
	dat->best_seed_idx = malloc(opt->K * sizeof *dat->best_seed_idx);
	dat->seeds = malloc(opt->K * sizeof *dat->seeds);
	dat->best_modes = malloc(opt->K * sizeof *dat->best_modes);

	if (dat->best_seed_idx == NULL || dat->seeds == NULL
		|| dat->best_modes == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::seeds");
	
	if (opt->filter_haplotypes) {
		dat->ini_seeds = malloc(opt->K * sizeof *dat->ini_seeds);
		dat->ini_seed_idx = malloc(opt->K * sizeof *dat->ini_seed_idx);
		if (!dat->ini_seeds || !dat->ini_seed_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"data::ini_*");
	}

	data_t *tmp1 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
	data_t *tmp2 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
	data_t *tmp3 = NULL;

	if (!tmp1 || !tmp2)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::seeds");

	if (opt->filter_haplotypes) {
		tmp3 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
		if (!tmp3)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"data::ini_seeds");
	}

	for (unsigned int i = 0; i < opt->K; i++) {
		dat->seeds[i] = tmp1;
		dat->best_modes[i] = tmp2;
		tmp1 += dat->n_coordinates;
		tmp2 += dat->n_coordinates;
		if (opt->filter_haplotypes) {
			dat->ini_seeds[i] = tmp3;
			tmp3 += dat->n_coordinates;
		}
	}

	return NO_ERROR;
} /* make_seeds */

/**
 * Setup data object.
 */

int make_hap_data(data **dat)
{
	*dat = malloc(sizeof **dat);

	if (*dat == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");

	(*dat)->fdata = NULL;
	(*dat)->dmat = NULL;
	(*dat)->qmat = NULL;
	(*dat)->n_quality = 0;
	(*dat)->data = NULL;
	
	(*dat)->read_idx = NULL;
	(*dat)->max_read_length = 0;
	(*dat)->min_read_length = 0;
	(*dat)->lengths = NULL;

	(*dat)->true_cluster_id = NULL;
	(*dat)->true_cluster_size = NULL;

	(*dat)->categories = NULL;
	(*dat)->n_observations = 0;
	(*dat)->n_coordinates = 0;
	(*dat)->n_categories = NULL;
	(*dat)->max_n_categories = 0;
	(*dat)->n_categories = NULL;
	(*dat)->tot_n_categories = 0;

	(*dat)->seeds = NULL;
	(*dat)->ini_seeds = NULL;
	(*dat)->seed_idx = NULL;
	(*dat)->best_seed_idx = NULL;
	(*dat)->ini_seed_idx = NULL;

	(*dat)->total = 0;
	(*dat)->cluster_id = NULL;
	(*dat)->obsn_idx = NULL;
	(*dat)->best_obsn_idx = NULL;
	(*dat)->best_cluster_id = NULL;
	(*dat)->best_modes = NULL;
	(*dat)->criterion = NULL;
	(*dat)->best_criterion = NULL;
	(*dat)->cluster_size = NULL;
	(*dat)->best_cluster_size = NULL;
	(*dat)->best_seed_idx = NULL;

	(*dat)->best_total = -INFINITY;
	(*dat)->best_rand = -INFINITY;
	(*dat)->use_ini = 0;
	(*dat)->n_init = 0;
	(*dat)->iter = 0;
	(*dat)->max_iter = 1;
	(*dat)->seconds = 0.;
	(*dat)->uncounted_seconds = 0.;
	(*dat)->first_cost = 1;
	(*dat)->worst_cost = 1;
	(*dat)->avg_cost = (*dat)->sd_cost = 0;
	(*dat)->avg_iter = (*dat)->sd_iter = 0;
	(*dat)->avg_time = (*dat)->sd_time = 0;
	(*dat)->ctime = 0;
	(*dat)->ntimes = 0;
	(*dat)->avg_ar = (*dat)->sd_ar = 0;
	(*dat)->avg_mi = (*dat)->sd_mi = 0;
	(*dat)->avg_vi = (*dat)->sd_vi = 0;

	(*dat)->seq_count = NULL;
	(*dat)->hash_length = 0;
	(*dat)->hash_leng_exk = 0;
	(*dat)->hash_idx_exk= NULL;
	(*dat)->hash_idx_exk_arr = NULL;
	(*dat)->llk_seq = NULL;
	(*dat)->sel_uniq_seq_count = NULL;
	(*dat)->min_index = NULL;

	(*dat)->uniq_seq_idx = NULL;
	(*dat)->uniq_seq_count = NULL;
	(*dat)->reads_uniq_id = NULL;
	(*dat)->exp_no_errs = NULL;
	(*dat)->mean_exp_no_errs = NULL;
	(*dat)->idx_arr = NULL;
	(*dat)->sel_uniq_seq_idx = NULL;
	(*dat)->sel_uniq_seq_count = NULL;
	(*dat)->sel_uniq_seq_idx = NULL;
	(*dat)->n_selected_sequences = 0;

	return NO_ERROR;
}/* make_data */

void free_hap_data(data *dat)
{
	if (!dat) return;

	if (dat->fdata) free_fastq(dat->fdata);
	if (dat->dmat) free(dat->dmat);
	if (dat->qmat) free(dat->qmat);
	if (dat->read_idx) free(dat->read_idx);
	if (dat->lengths) free(dat->lengths);
	if (dat->true_cluster_id) free(dat->true_cluster_id);
	if (dat->true_cluster_size) free(dat->true_cluster_size);
	
	if (dat->seeds) {
		free(dat->seeds[0]);
		free(dat->seeds);
	}
	if (dat->ini_seeds) {
		free(dat->ini_seeds[0]);
		free(dat->ini_seeds);
	}
	if (dat->best_modes) {
		free(dat->best_modes[0]);
		free(dat->best_modes);
	}

	if (dat->seed_idx) free(dat->seed_idx);
	if (dat->best_seed_idx) free(dat->best_seed_idx);
	if (dat->cluster_id) free(dat->cluster_id);
	if (dat->best_cluster_id) free(dat->best_cluster_id);
	if (dat->criterion) free(dat->criterion);
	if (dat->best_criterion) free(dat->best_criterion);
	if (dat->cluster_size) free(dat->cluster_size);
	if (dat->best_cluster_size) free(dat->best_cluster_size);
	if (dat->n_categories) free(dat->n_categories);
	if (dat->categories) free(dat->categories);

	if (dat->seq_count) delete_all(dat->seq_count);
	if (dat->hash_idx_exk) free(dat->hash_idx_exk);
	if (dat->hash_idx_exk_arr) free(dat->hash_idx_exk_arr);
	if (dat->llk_seq) free(dat->llk_seq);
	if (dat->min_index) free(dat->min_index);

	if (dat->uniq_seq_count) free(dat->uniq_seq_count);
	if (dat->uniq_seq_idx) free(dat->uniq_seq_idx);
	if (dat->reads_uniq_id) free(dat->reads_uniq_id);
	if (dat->mean_exp_no_errs) free(dat->mean_exp_no_errs);
	if (dat->exp_no_errs) free(dat->exp_no_errs);
	if (dat->sel_uniq_seq_idx) free(dat->sel_uniq_seq_idx);
	if (dat->sel_idx_arr) free(dat->sel_idx_arr);
	if (dat->sel_uniq_seq_count) free(dat->sel_uniq_seq_count);

	free(dat);

} /* free_data */

/**
 * Write solution
 *
 * @param dat		data object pointer
 * @param opt		options object pointer
 * @param in_fps	solution file pointer
 */
int write_solution(data *dat, options *opt, FILE **in_fps)
{
	/* report best solution to stdout and outfile */

	*in_fps = fopen(opt->soln_file, "w");
	if (*in_fps == NULL)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->soln_file);
	FILE *fps = *in_fps;

	if (fps && opt->K > 0) {
		fprintf(fps, "Seed: %lu\n", opt->seed);
		fprintf(fps, "Algorithm: %s\n",
			khaplotype_algorithm(opt->khaplotype_algorithm));
		fprintf(fps, "Initialization: %s%s\n",
			khaplotype_init_method(opt->init_method),
			opt->shuffle ? " with shuffling." : ".");
		fprintf(fps, "Run for %u initializations\n", opt->n_init);
		fprintf(fps, "Input file: %s\n", opt->datafile);
		fprintf(fps, "Number of clusters to estimate: %u\n", opt->K);
		fprintf(fps, "Number of observations: %u\n",
			dat->n_observations);
		fprintf(fps, "Number of coordinates: %u\n", dat->n_coordinates);
		if (opt->run_with_quals == 0) {
			fprintf(fps, "Number of categories: ");
			fprint_data_ts(fps, dat->n_categories, dat->n_coordinates,
					(int)(log10(dat->max_n_categories -
							(dat->max_n_categories>1))+1), 1);
		}
		fprintf(fps, "Total categories: %hhu\n", dat->tot_n_categories);
		fprintf(fps, "Running time: %0.6lf\n", dat->seconds);
	}

	if (fps) {
		fprintf(fps, "Best optimized criterion: %.0f\n",
			dat->best_total);
		fprint_doubles(fps, dat->best_criterion, opt->K, 0, 1);
	}

	if (fps && opt->K > 0) {
		fprintf(fps, "Best cluster sizes:");
		fprint_uints(fps, dat->best_cluster_size, opt->K, 0, 1);
	}

	if (fps && opt->K > 0) {
		fprintf(fps, "Best solution cluster assignments:\n");
		fprint_uints(fps, dat->best_cluster_id,
						dat->n_observations, 0, 1);
		fprintf(fps, "\n");
	}

	if (fps) {
		fprintf(fps, "Best modes:\n");
		for (unsigned int k = 0; k < opt->K; ++k)
			fprint_data_ts(fps, dat->best_modes[k],
						dat->n_coordinates, 0, 1);
	}

	fclose(fps);
	return NO_ERROR;
} /* write_solution */

/**
 * Shuffle data input order.  Since the k-modes algorithms take the observations
 * in input order for initialization, and many of the real datasets are ordered
 * by class, it is important to randomize this aspect of the data.
 *
 * See https://stackoverflow.com/questions/3343797/is-this-c-implementation-of-fisher-yates-shuffle-correct
 * and https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 *
 * @param dat	data pointer
 * @return	error status
 */
int shuffle_data_hap(data *dat, options *opt)
{
	if (dat->n_observations > RAND_MAX)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Data set too big!");
	long n = (long) dat->n_observations;
	long i, j, lim;
	data_t *dptr, *qptr;

	for (i = n - 1; i > 0; --i) {
		lim = RAND_MAX - RAND_MAX % (i + 1);
		do {
			j = rand();
		} while (j >= lim);
		j = j % (i + 1);

		dptr = dat->dmat[j];
		dat->dmat[j] = dat->dmat[i];
		dat->dmat[i] = dptr;

		if (opt->run_with_quals) {
			qptr = dat->qmat[j];
			dat->qmat[j] = dat->qmat[i];
			dat->qmat[i] = qptr;
		}

		if (dat->obsn_idx) {
			unsigned int ui = dat->obsn_idx[j];
			dat->obsn_idx[j] = dat->obsn_idx[i];
			dat->obsn_idx[i] = ui;
		}
	}

	return NO_ERROR;
} /* shuffle_data */

/**
 * Build hash of reads.
 *
 * @param dat	pointer to data object
 * @return	error status
 */
int build_hash(data *dat, options *opt)
{
	int err = NO_ERROR;
	/* build hash table */
	dat->hash_length = 0;
	for (size_t i = 0; i < dat->n_observations; ++i) {
		//		fprintf(stderr, "%zu: %s\n", i, display_sequence(dat->dmat[i], dat->lengths[i], dat->fdata->read_encoding));
		dat->hash_length += add_sequence(&dat->seq_count, dat->dmat[i],
						 dat->lengths[i], i);
	}

	/* store index of reads for all unique sequences */
	for (size_t i = 0; i < dat->n_observations; ++i)
		if ((err = add_seq_idx(dat->seq_count, dat->dmat[i],
							dat->lengths[i], i)))
			return err;

	/* sort hash table by count */
	sort_by_count(&dat->seq_count);

	/* allocate for hash with abundance more than k */
	dat->hash_leng_exk = count_sequences(dat->seq_count, opt->hash_k);
	MAKE_1ARRAY(dat->hash_idx_exk, dat->hash_leng_exk);

	MAKE_1ARRAY(dat->seq_abundance, dat->hash_leng_exk);
	MAKE_1ARRAY(dat->min_index, opt->K);

	/* store the index of unique sequence with abundance exceed k in the hash */
	store_idx_exk(dat->seq_count, dat->hash_idx_exk, dat->hash_length,
		opt->hash_k, &dat->hash_idx_exk_arr, dat->seq_abundance);

	return err;
} /* build_hash */

/**
 * Find haplotypes in a hash table with mean expected error exceeding cut-off
 * value.  This function could be adapted to filter haplotypes by other methods.
 *
 * @param dat	data sturcture
 * @param opt	option
 * @return	error status
 */
int
filter_haplotypes (data *dat, options *opt)
{

	dat->exp_no_errs = malloc(dat->n_observations
					* sizeof *dat->exp_no_errs);
	dat->mean_exp_no_errs = calloc(dat->hash_length,
					sizeof *dat->mean_exp_no_errs);

	if (!dat->exp_no_errs || !dat->mean_exp_no_errs)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data.exp_no_errs");

	/* allocate room for index and count array of unique sequences */
	dat->uniq_seq_idx = calloc(dat->hash_length, sizeof *dat->uniq_seq_idx);
	dat->uniq_seq_count = calloc(dat->hash_length,
						sizeof *dat->uniq_seq_count);
	dat->reads_uniq_id = calloc(dat->n_observations,
						sizeof *dat->reads_uniq_id);

	if (!dat->uniq_seq_idx || !dat->uniq_seq_count || !dat->reads_uniq_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"initializer.uniq_seq");

	dat->sel_uniq_seq_idx = malloc(dat->hash_length
						* sizeof *dat->sel_uniq_seq_idx);
	dat->sel_uniq_seq_count = malloc(dat->hash_length
						* sizeof *dat->sel_uniq_seq_count);
	dat->sel_idx_arr = malloc(dat->hash_length
						* sizeof *dat->sel_idx_arr);

	/* sort hash table in an order of decreasing abundance */
	if (store_index(dat->seq_count, dat->reads_uniq_id, dat->uniq_seq_idx,
		&dat->idx_arr, dat->hash_length, dat->n_observations)
		|| store_count(dat->seq_count, dat->uniq_seq_count,
		dat->hash_length))
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
				"store_index() or store_count()");

	for (unsigned int i = 0; i < dat->n_observations; ++i) {
		dat->exp_no_errs[i] = 0;
		for (unsigned int j = 0; j < dat->max_read_length; ++j)
			dat->exp_no_errs[i]
				+= error_prob(dat->fdata, dat->qmat[i][j]);
		dat->mean_exp_no_errs[dat->reads_uniq_id[i]]
							+= dat->exp_no_errs[i];
	}

	for (unsigned int i = 0; i < dat->hash_length; ++i) {
		dat->mean_exp_no_errs[i] /= dat->uniq_seq_count[i];
		if (dat->mean_exp_no_errs[i] < opt->cutoff_mean_exp_err) {
			dat->sel_uniq_seq_idx[dat->n_selected_sequences]
							= dat->uniq_seq_idx[i];
			dat->sel_uniq_seq_count[dat->n_selected_sequences]
							= dat->uniq_seq_count[i];
			dat->sel_idx_arr[dat->n_selected_sequences]
							= dat->idx_arr[i];
			dat->n_selected_sequences++;
		}
	}

	return NO_ERROR;
}/* filter_haplotypes */

/**
 * Pre-compute log probabilities of no error or error given quality scores.
 *
 * @param d		void pointer to data
 * @param prob_t	(ln(1-pij))
 * @param prob_f	(ln (pij/3))
 */
void
compute_pij (void *d, double **prob_t, double **prob_f)
{
	data *dp = (data *)d;
	unsigned int i, j;
	double prob;
	int flag_q = dp->qmat != NULL;
	
	/* To see if include the quality scores */
	if (flag_q) {
		/* Store the prob for each read if have quality info */
		for (i = 0; i < dp->n_observations; ++i)
			for (j = 0; j < dp->lengths[i]; ++j) {
				prob = error_prob(dp->fdata, dp->qmat[i][j]);
				prob_t[i][j] = log(1 - prob);
				prob_f[i][j] = log(prob) - LOG3;
			}
	} else {
		/* Store the prob for each read if no quality info */
		for (i = 0; i < dp->n_observations; ++i)
			for (j = 0; j < dp->n_coordinates; ++j) {
				prob_t[i][j] = 1;
				prob_f[i][j] = 0;
			}
	}
}

int read_hap_data(data *dat, options *opt)
{
	int err = NO_ERROR;
	unsigned int j = 0, i;
	char c, pc = 0;
	FILE *fp = fopen(opt->datafile, "r");
	data_t *dptr;
	
	if (fp == NULL)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->datafile);
	
	dat->n_coordinates = 1;
	dat->n_observations = 0;
	
	/* count number of columns, assuming space-separated fields,
	 * allowing possibility of no terminal newline */
	while (!feof(fp)) {
		c = fgetc(fp);
		/* new line starting */
		if (!feof(fp) && c != '\n' && (!pc || pc == '\n'))
			dat->n_observations++;
		/* new column starting */
		if (dat->n_observations == 1 && pc && c == ' ' && pc != ' ')
			dat->n_coordinates++;
		pc = c;
	}
	
	fprintf(stderr, "quiet = %d\n", opt->quiet);
	debug_msg(MINIMAL, opt->quiet, "Data %u x %u\n", dat->n_observations,
		  dat->n_coordinates);
	
	if (opt->true_column < dat->n_coordinates) {
		--dat->n_coordinates;
		opt->true_cluster = malloc(dat->n_observations
					   * sizeof *opt->true_cluster);
		if (!opt->true_cluster) {
			err = mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				       "options::true_cluster");
			goto ABORT_READ_DATA;
		}
	}
	
	dat->n_categories = calloc(dat->n_coordinates,
				   sizeof *dat->n_categories);
	dat->data = malloc(dat->n_coordinates * dat->n_observations
			   * sizeof *dat->data);
	if (dat->data == NULL || dat->n_categories == NULL) {
		err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::data");
		goto ABORT_READ_DATA;
	}
	
	dat->cluster_id = malloc(dat->n_observations * sizeof *dat->cluster_id);
	dat->best_cluster_id = malloc(dat->n_observations
				      * sizeof *dat->best_cluster_id);
	
	if (dat->cluster_id == NULL || dat->best_cluster_id == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::cluster_id");
	
	rewind(fp);
	
	/* read in data */
	/* assume categories are 0, 1, 2, ... without skips */
	dptr = dat->data;
	if (opt->true_column > dat->n_coordinates) {
		unsigned int nc = dat->n_coordinates;
		while (fscanf(fp, "%" AMPLICLUST_SCNu_data_t, dptr) == 1) {
			if (opt->subtract_one) -- (*dptr);
			if (dat->n_categories[j] < *dptr + 1u)
				dat->n_categories[j] = *dptr + 1u;
			++ dptr;
			j = (j + 1u) % nc;
		}
	} else {
		unsigned int nc = dat->n_coordinates + 1u;
		i = 0;
		while (fscanf(fp, "%" AMPLICLUST_SCNu_data_t, dptr) == 1) {
			if (opt->subtract_one) -- (*dptr);
			if (j == opt->true_column) {
				opt->true_cluster[i ++] = *dptr;
				if (*dptr + 1u > opt->true_K)
					opt->true_K = *dptr + 1u;
			} else {
				if (*dptr + 1u > dat->n_categories[j
								   - (j > opt->true_column)])
					dat->n_categories[j - (j
							       > opt->true_column)]
					= *dptr + 1u;
				++ dptr;
			}
			j = (j + 1u) % nc;
		}
	}
	fclose(fp);
	fp = NULL;
	
	/* fix assumption of 0, 1, 2, ... categories */
	/* SLOW: should write an option to bypass it! */
	for (j = 0; j < dat->n_coordinates; ++j) {
		unsigned int *present = calloc(dat->n_categories[j], sizeof *present);
		for (i = 0; i < dat->n_observations; ++i)
			++present[dat->data[dat->n_coordinates*i + j]];
		unsigned int cnt = 0;
		for (i = 0; i < dat->n_categories[j]; ++i) {
			if (present[i]) ++cnt;
			present[i] = (i ? present[i-1] : 0) + (present[i] == 0);
		}
		if (cnt == dat->n_categories[j]) {
			free(present);
			continue;
		}
		mmessage(WARNING_MSG, NO_ERROR, "Coordinate %u uses only %u "
			 "categories, but has %u categories.  To get correct "
			 "category counts, use double arguments to -o and maybe "
			 "-m, and then run with corrected files.\n", j,
			 cnt, dat->n_categories[j]);
		dat->n_categories[j] = cnt;
		for (i = 0; i < dat->n_observations; ++i)
			dat->data[dat->n_coordinates*i + j] -= present[dat->data[dat->n_coordinates*i + j]];
		free(present);
	}
	
	for (j = 0; j < dat->n_coordinates; ++j) {
		dat->tot_n_categories += dat->n_categories[j];
		if (dat->n_categories[j] > dat->max_n_categories)
			dat->max_n_categories = dat->n_categories[j];
	}
	
ABORT_READ_DATA:
	if (fp) fclose(fp);
	
	return err;
} /* read_data */

int finish_make_hap_data(data *dat, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	
	dat->dmat = malloc(dat->n_observations * sizeof *dat->dmat);
	data_t *rptr = dat->data;
	for (unsigned int i = 0; i < dat->n_observations; i++) {
		dat->dmat[i] = rptr;
		rptr += dat->n_coordinates;
	}
	if (fxn_debug)
		mmessage(DEBUG_MSG, NO_ERROR, "Allocated %dx%d data matrix\n",
			 dat->n_observations, dat->n_coordinates);
	
	/* Find the number of category and count it */
	data_t n_unique;
	dat->categories = malloc(dat->tot_n_categories * sizeof *dat->categories);
	if (!dat->categories)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.cat");\
	
	dat->categories = find_unique(dat->data, dat->n_observations * dat->n_coordinates, &n_unique);
	dat->tot_n_categories = n_unique;
	
	/* Allocate for writing results */
	dat->cluster_id = malloc(dat->n_observations * sizeof *dat->cluster_id);
	dat->best_cluster_id = malloc(dat->n_observations * sizeof *dat->best_cluster_id);
	
	if (dat->cluster_id == NULL || dat->best_cluster_id == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::cluster_id");
	
	//    MAKE_2ARRAY(dat->best_modes, opt->K, dat->max_read_length);
	
	dat->cluster_size = malloc(opt->K * sizeof *dat->cluster_size);
	dat->best_cluster_size = malloc(opt->K * sizeof *dat->best_cluster_size);
	
	if (dat->cluster_size == NULL || dat->best_cluster_size == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::cluster_size");
	
	dat->criterion = malloc(opt->K * sizeof *dat->criterion);
	dat->best_criterion = malloc(opt->K * sizeof *dat->best_criterion);
	
	if (dat->criterion == NULL || dat->best_criterion == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::criterion");
	
	if (opt->shuffle) {
		dat->obsn_idx = malloc(dat->n_observations
				       * sizeof *dat->obsn_idx);
		dat->best_obsn_idx = malloc(dat->n_observations
					    * sizeof *dat->obsn_idx);
		if (!dat->obsn_idx || !dat->best_obsn_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"data:obsn_idx");
		for (unsigned int i = 0; i < dat->n_observations; ++i)
			dat->obsn_idx[i] = i;
	}
	
	if (!opt->K)
		return NO_ERROR;
	
	if (!opt->seed_idx) {
		dat->seed_idx = malloc(opt->K * sizeof *dat->seed_idx);
		if (!dat->seed_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"data::seeds_idx");
	} else
		dat->seed_idx = opt->seed_idx;
	
	dat->best_seed_idx = malloc(opt->K * sizeof *dat->best_seed_idx);
	dat->seeds = malloc(opt->K * sizeof *dat->seeds);
	dat->best_modes = malloc(opt->K * sizeof *dat->best_modes);
	
	if (dat->best_seed_idx == NULL || dat->seeds == NULL
	    || dat->best_modes == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::seeds");
	
	data_t *tmp1 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
	data_t *tmp2 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
	if (!tmp1 || !tmp2)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::seeds");
	for (unsigned int i = 0; i < opt->K; i++) {
		dat->seeds[i] = tmp1;
		dat->best_modes[i] = tmp2;
		tmp1 += dat->n_coordinates;
		tmp2 += dat->n_coordinates;
	}
	return NO_ERROR;
}

int make_res(outres **out, data *dat, options *opt) {
	*out = malloc(sizeof **out);
	if (!*out)
		return(mmessage(ERROR_MSG, MEMORY_ALLOCATION, "out_data"));
	(*out)->n_observations = 0;
	(*out)->n_coordinates = 0;
	(*out)->best_total = 0;
	(*out)->best_criterion = NULL;
	(*out)->best_cluster_id = NULL;
	(*out)->best_cluster_size = NULL;
	(*out)->best_modes = NULL;
	
	MAKE_1ARRAY((*out)->best_criterion, opt->K);
	MAKE_1ARRAY((*out)->best_cluster_size, opt->K);
	MAKE_1ARRAY((*out)->best_cluster_id, dat->n_observations);
	MAKE_1ARRAY((*out)->best_modes, opt->K * dat->n_coordinates);
	
	(*out)->n_coordinates = dat->n_coordinates;
	(*out)->n_observations = dat->n_observations;
	(*out)->best_total = dat->best_total;
	COPY_1ARRAY((*out)->best_criterion, dat->best_criterion, opt->K);
	COPY_1ARRAY((*out)->best_cluster_size, dat->best_cluster_size, opt->K);
	COPY_1ARRAY((*out)->best_cluster_id, dat->best_cluster_id, dat->n_observations);
	for (unsigned int i = 0; i < opt->K; ++i)
		for (unsigned int j = 0; j < dat->n_coordinates; ++j)
			(*out)->best_modes[i * dat->n_coordinates + j] = dat->best_modes[i][j];
	
	return NO_ERROR;
}

void free_res(outres *out)
{
	if (out->best_cluster_id) free(out->best_cluster_id);
	if (out->best_criterion) free(out->best_criterion);
	if (out->best_cluster_size) free(out->best_cluster_size);
	if (out->best_modes) free(out->best_modes);
}
