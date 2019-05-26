/**
 * @file run_khaplotype.c
 * @author Yudi Zhang
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>

#include "khaplotype.h"
#include "khaplotype_data.h"
#include "khaplotype_option.h"
#include "khaplotype_model.h"
#include "khaplotype_initialize.h"
#include "cluster_lloyds.h"
#include "array.h"
#include "cluster.h"
#include "fastq.h"
#include "cmdline.h"
#include "timing_mach.h"
#include "error.h"

static inline void stash_state(data *dat, options *opt);
int khaplotype(options *opt, outres **out);

#ifndef STANDALONE
#include <Rinternals.h>
#define PRINTF(str, ...) Rprintf((str), __VA_ARGS__)
#define EPRINTF(str, ...) REprintf((str), __VA_ARGS__)
#else
#define PRINTF(str, ...) fprintf(stdout, (str), __VA_ARGS__)
#define EPRINTF(str, ...) fprintf(stderr, (str), __VA_ARGS__)
#endif


#ifdef STANDALONE
int main(int argc, const char **argv)
{
	int err = NO_ERROR;
	options *opt = NULL;		/* command-line options */
	outres *out = NULL;		/* data object */

	/* parse command line */
	if ((err = make_opt(&opt)))
		goto CLEAR_AND_EXIT;
	
	if ((err = parse_opt(opt, argc, argv)))
		goto CLEAR_AND_EXIT;
	
	khaplotype(opt, &out);

//	printf("%d\n", out->n_observations);
//	printf("%d\n", out->n_coordinates);
//	printf("%lf\n", out->best_total);
//	for (unsigned int k = 0; k < opt->K; ++k)
//		fprint_data_ts(stderr, out->best_modes[k],
//			       out->n_coordinates, 0, 1);
//	fprint_uints(stderr, out->best_cluster_id,
//		     out->n_observations, 0, 1);
	
CLEAR_AND_EXIT:

	if (opt)
		free_opt(opt);
	if (out)
		free_res(out);
	
	return(err);
} /* main */

#else
SEXP r_khaplotype (SEXP K_r,
		   SEXP datafile_r,
		   SEXP n_init_r,
		   SEXP khaplotype_algorithm_r,
		   SEXP seed_r,
//		   SEXP sim_info_file_r,
		   SEXP shuffle_r,
		   SEXP run_with_quals_r
		   )
{
	int err = NO_ERROR;
	options *opt = NULL;		/* command-line options */
	outres *out = NULL;		/* data object */
	SEXP r_list;
	/* parse command line */
	if ((err = make_opt(&opt)))
		return R_NilValue;

	opt->K = *INTEGER(K_r);
	opt->n_init = *INTEGER(n_init_r);
	opt->shuffle = *INTEGER(shuffle_r);
	opt->run_with_quals = *INTEGER(run_with_quals_r);
	opt->seed = *INTEGER(seed_r);
	srand(opt->seed);
	opt->khaplotype_algorithm = *INTEGER(khaplotype_algorithm_r);
	opt->datafile = CHAR(STRING_ELT(datafile_r, 0));
//	opt->sim_info_file = isNull(sim_info_file_r) ? NULL : CHAR(STRING_ELT(sim_info_file_r, 0));

	khaplotype(opt, &out);
	
	SEXP r_best_cluster_size = PROTECT(allocVector(INTSXP, opt->K));
	memcpy(INTEGER(r_best_cluster_size), out->best_cluster_size, opt->K * sizeof *out->best_cluster_size);
	
	SEXP r_best_cluster_id = PROTECT(allocVector(INTSXP, out->n_observations));
	memcpy(INTEGER(r_best_cluster_id), out->best_cluster_id, out->n_observations * sizeof *out->best_cluster_id);
	
	SEXP r_best_modes = PROTECT(allocVector(INTSXP, opt->K * out->n_coordinates));
	memcpy(INTEGER(r_best_modes), out->best_modes, out->n_coordinates * opt->K * sizeof *out->best_modes);
	
	SEXP r_data = PROTECT(allocVector(INTSXP, out->n_observations * out->n_coordinates));
	memcpy(INTEGER(r_data), out->data, out->n_coordinates * out->n_observations * sizeof *out->data);

	SEXP r_best_criterion = PROTECT(allocVector(REALSXP, opt->K));
	memcpy(REAL(r_best_criterion), out->best_criterion, opt->K * sizeof *out->best_criterion);
	
	PROTECT(r_list = allocVector(VECSXP, 5));
	SET_VECTOR_ELT(r_list, 0, r_best_cluster_size);
	SET_VECTOR_ELT(r_list, 1, r_best_criterion);
	SET_VECTOR_ELT(r_list, 2, r_best_cluster_id);
	SET_VECTOR_ELT(r_list, 3, r_best_modes);
	SET_VECTOR_ELT(r_list, 4, r_data);

	if (opt)
		free_opt(opt);
	if (out)
		free_res(out);
	UNPROTECT(6);
	return r_list;
}

SEXP r_read_fastq(SEXP datafile_r)
{
	int err = NO_ERROR;
	fastq_options *fqo = NULL;
	fastq_data *fdata = NULL;
	char const *datafile = NULL;
	SEXP r_list;
	
	
	datafile = CHAR(STRING_ELT(datafile_r, 0));
	
	if ((err = make_fastq_options(&fqo)))
		return R_NilValue;
	fqo->read_encoding = XY_ENCODING;
	read_fastq(datafile, &fdata, fqo);
	
	SEXP r_data_dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(r_data_dim)[0] = fdata->n_reads;
	INTEGER(r_data_dim)[1] = fdata->n_max_length;
	
	SEXP r_reads = PROTECT(allocVector(INTSXP, fdata->n_reads * fdata->n_max_length));
	SEXP r_quals = PROTECT(allocVector(INTSXP, fdata->n_reads * fdata->n_max_length));
	
	for (unsigned int i = 0; i < fdata->n_reads * fdata->n_max_length; ++i) {
		INTEGER(r_reads)[i] = fdata->reads[i];
		INTEGER(r_quals)[i] = fdata->quals[i];
	}
	
	PROTECT(r_list = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(r_list, 0, r_reads);
	SET_VECTOR_ELT(r_list, 1, r_quals);
	SET_VECTOR_ELT(r_list, 2, r_data_dim);
	UNPROTECT(4);
	
	if (fdata)
		free_fastq(fdata);
	
	return r_list;
}
#endif

int khaplotype(options *opt, outres **out)
{
	int err = NO_ERROR;
	data *dat = NULL;		/* data object */
	auxvar *aux = NULL;		/* auxiliary variables structure */
	model *mod = NULL;
	fastq_options *fqo = NULL;	/* fastq file options */
	double ari = 0, best_ari = -1;
	TIME_STRUCT start, each_run;		/* timing */
	
	/* make data and auxiliary variables object */
	if ((err = make_hap_data(&dat)))
		goto CLEAR_AND_EXIT;
	if ((err = make_aux_var(&aux)))
		goto CLEAR_AND_EXIT;
	
	if (opt->run_with_quals) {
		
		if ((err = make_fastq_options(&fqo)))
			goto CLEAR_AND_EXIT;
		
		/* encode nucleotides in 2-bits: error if ambiguous bases */
		fqo->read_encoding = XY_ENCODING;
		
		/* read sequence data */
		if (opt->datafile && (err = read_fastq(opt->datafile,
						       &dat->fdata, fqo)))
			goto CLEAR_AND_EXIT;
		
		/* with data now loaded, can polish off data object */
		if ((err = sync_state_khaplotype_data(dat, opt)))
			goto CLEAR_AND_EXIT;
		if ((err = sync_state_aux_var(dat, aux, opt)))
			goto CLEAR_AND_EXIT;
		
		/* build hash table to find unique sequence */
		if ((err = build_hash(dat, opt)))
			goto CLEAR_AND_EXIT;
		
		/* optionally remove some haplotypes from seeding */
		if (opt->filter_haplotypes)
			if ((err = filter_haplotypes(dat, opt)))
				return err;
	}
	else {
		/* run on normal categorical data */
		if ((err = read_hap_data(dat, opt)))
			goto CLEAR_AND_EXIT;
		if ((err = finish_make_hap_data(dat, opt)))
			goto CLEAR_AND_EXIT;
		if ((err = sync_state_aux_var(dat, aux, opt)))
			goto CLEAR_AND_EXIT;
	}
	
	/* prepare seeds */
	if ((err = make_seeds(dat, opt)))
		goto CLEAR_AND_EXIT;
	
	/* model object mod for storing best result under different K */
	if (opt->init_method == KHAPLOTYPE_INIT_FILTER_DATA)
		if ((err = make_model(&mod, dat)))
			return err;
	
#ifdef STANDALONE
	FILE *fps = NULL;
#endif
	
	/* Record time */
	MARK_TIME(&start);
	
	for (unsigned int i = 0; i < opt->n_init; ++i) {
		
		/* shuffle data to avoid input order artifacts */
		if (opt->shuffle && opt->K > 1 && (err = shuffle_data_hap(dat, opt)))
			goto CLEAR_AND_EXIT;
		
		/* initialization */
		if ((err = initialize(dat, aux, opt, mod)))
			goto CLEAR_AND_EXIT;
		
		MARK_TIME(&each_run);
		/* Run one of the algorithms */
		if (opt->khaplotype_algorithm == FASTQ_LLOYDS) {
			dat->total = cluster_lloyds(
						    dat, aux,
						    dat->use_ini ? dat->ini_seeds : dat->seeds,
						    dat->cluster_size, dat->cluster_id,
						    dat->n_observations, dat->n_coordinates, opt->K,
						    opt->n_max_iter, dat->criterion,
						    &err, &dat->iter,
						    fastq_lloyds_step1, fastq_lloyds_step2,
						    compute_criterion_hap);
			
		} else if (opt->khaplotype_algorithm == FASTQ_LLOYDS_EFFICIENT) {
			dat->total = cluster_lloyds_efficient(
							      dat, aux,
							      dat->use_ini ? dat->ini_seeds : dat->seeds,
							      dat->cluster_size, dat->cluster_id,
							      dat->n_observations, dat->n_coordinates, opt->K,
							      opt->n_max_iter, dat->criterion,
							      &err, &dat->iter,
							      fastq_lloyds_efficient_step1,
							      fastq_lloyds_efficient_step2,
							      compute_cost);
			
		} else if (opt->khaplotype_algorithm == FASTQ_MACQUEEN) {
			dat->total = cluster_macqueen(
						      dat, aux,
						      dat->n_observations, dat->n_coordinates,
						      dat->use_ini ? dat->ini_seeds : dat->seeds,
						      opt->K,
						      dat->cluster_id, dat->cluster_size,
						      opt->n_max_iter, dat->criterion,
						      &err, &dat->iter,
						      fastq_macqueen_ini, fastq_macqueen_iter,
						      compute_cost);
		} else if (opt->khaplotype_algorithm == FASTQ_HW) {
			dat->total = cluster_hw(
						dat, aux,
						dat->use_ini ? dat->ini_seeds : dat->seeds,
						dat->cluster_size, dat->cluster_id,
						dat->n_observations, dat->n_coordinates, opt->K,
						opt->n_max_iter, dat->criterion,
						&err, &dat->iter,
						hw_init_default, hw_iter_default);
			
		} else if (opt->khaplotype_algorithm == FASTQ_HW_EFFICIENT) {
			dat->total = cluster_hw_fast(
						     dat, aux,
						     dat->use_ini ? dat->ini_seeds : dat->seeds,
						     dat->cluster_size, dat->cluster_id,
						     dat->n_observations, dat->n_coordinates, opt->K,
						     opt->n_max_iter, dat->criterion,
						     &err, &dat->iter,
						     hw_fast_init_default, hw_fast_iter_default);
		} else {
#ifndef STANDALONE
			warning("[WARNING] Invalid algorithm.\n");
#else
			mmessage(ERROR_MSG, INVALID_USER_INPUT,
				 "Invalid algorithm.\n");
#endif
			goto CLEAR_AND_EXIT;
		}
		dat->uncounted_seconds = ELAP_TIME(&each_run);
		
		if (err) {
			mmessage(ERROR_MSG, CUSTOM_ERROR, "%s\n",
				 khaplotype_error(err));
			if (err != KHAPLOTYPE_EXCEED_ITER_WARNING)
				goto CLEAR_AND_EXIT;
		}
		
		/* record best solution */
		PRINTF("Time cost: %lf secs\n", dat->uncounted_seconds);
		PRINTF("Log likelihood in %dth initialization: %0.2lf (%u iterations: %u %u)\n", i + 1, dat->total, dat->iter, dat->cluster_size[0], dat->cluster_size[1]);
		if (opt->sim_cluster) {
			ari = cluster_index(dat->cluster_id,
					    opt->sim_cluster, dat->n_observations,
					    opt->K, opt->sim_K, ADJUSTED_RAND_INDEX);
			PRINTF("ARI: %f\n", ari);
		}
		if ((!err || err == KHAPLOTYPE_EXCEED_ITER_WARNING)
		    && dat->total > dat->best_total) {
			stash_state(dat, opt);
			if (opt->sim_cluster)
				best_ari = ari;
		} else if (err)
			PRINTF("[ERROR] %s (%d)\n",
				khaplotype_error(err), err);
	}
	dat->seconds = ELAP_TIME(&start);
	PRINTF("Time cost is: %lf secs\n", dat->seconds);
	PRINTF("Best optimum is: %lf\n", dat->best_total);
	if (opt->sim_cluster)
		PRINTF("ARI at optimum is: %f\n", best_ari);
	
	/* Write solution */
#ifndef STANDALONE
	make_res(out, dat, opt);
#else
	write_solution(dat, opt, &fps);
#endif
	
CLEAR_AND_EXIT:
	
	if (dat)
		free_hap_data(dat);
	if (aux)
		free_aux_var(aux);
	if (mod)
		free_model(mod);
	
	return(err);
}

/**
 * Stash current state in best_* slot.
 *
 * @param dat    data object pointer
 * @param opt    options object pointer
 */
static inline void stash_state(data *dat, options *opt)
{
	dat->best_total = dat->total;
	
	if (dat->use_ini)
		for (unsigned int k = 0; k < opt->K; ++k)
			memcpy(dat->best_modes[k], dat->ini_seeds[k],
			       dat->n_coordinates * sizeof **dat->best_modes);
	else
		for (unsigned int k = 0; k < opt->K; ++k)
			memcpy(dat->best_modes[k], dat->seeds[k], dat->n_coordinates * sizeof **dat->best_modes);
	
	COPY_1ARRAY(dat->best_cluster_id, dat->cluster_id, dat->n_observations);
	COPY_1ARRAY(dat->best_cluster_size, dat->cluster_size, opt->K);
	COPY_1ARRAY(dat->best_criterion, dat->criterion, opt->K);
	
	if (opt->shuffle)
		memcpy(dat->best_obsn_idx, dat->obsn_idx, dat->n_observations
		       * sizeof *dat->obsn_idx);
} /* stash_state */
