/**
 * @file haplotype_option.h
 * @author Yudi Zhang
 */

#ifndef KHAPLOTYPE_OPTION_H
#define KHAPLOTYPE_OPTION_H

#include <stdio.h>
#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>

#include "constants.h"

/** Varieties of algorithm for fastq data.
 */
enum {
	FASTQ_LLOYDS,
	FASTQ_LLOYDS_EFFICIENT,
	FASTQ_HW_EFFICIENT,
	FASTQ_MACQUEEN,
	FASTQ_HW
};

/**
 * Run options for k-haplotype.
 */
typedef struct _options options;

struct _options {
	/* model */
	unsigned int K;			/*<! number of clusters */
	double cut_k;			/*<! cut-off value when using BIC */
	
	/* data */
	int run_with_quals;		/*<! use quality scores */
	char const *datafile;		/*<! name of data file */
	int subtract_one;		/*<! subtract 1 from data categories */

	/* truth for simulated data */
	unsigned int true_column;	/*<! supervised data: column of truth */
	unsigned int *true_cluster;	/*<! true cluster assignments */
	unsigned int *true_cluster_size;/*<! true cluster sizes */
	unsigned int true_K;		/*<! true number of clusters */
	
	/* run conditions */
	unsigned int n_init;		/*<! no. random initializations */
	unsigned int n_inner_init;	/*<! no. random inner loop initializations */
	unsigned int n_max_iter;	/*<! max. number of iterations */
	unsigned long seed;		/*<! random number seed [srand()] */
	int khaplotype_algorithm;	/*<! run Lloyd's, Huang's or HW k-modes */
	int shuffle;			/*<! shuffle input order */
	
	/* initialization */
	int init_method;		/*<! initialization method to use */
	unsigned int *seed_idx;		/*<! user-provided seed indices */
	unsigned int n_sd_idx;		/*<! tmp: length of seed_idx */
	unsigned int n_seed_set;	/*<! number of seeds in seed set */
	unsigned int hash_k;		/*<! number of abundance in each unique sequence */
	unsigned int filter_haplotypes;	/*<! filter allowed haplotypes */
	double cutoff_mean_exp_err;	/*<! cut-off value for mean expected error */
	unsigned int run_with_hash;	/*<! use hash to pick the seeds */
	
	/* output */
	char const *soln_file;		/*<! solution file */
	int info;			/*<! level of information to output */
	int quiet;			/*<! be quiet */
	
	/* simulation */
	unsigned int sim_K;		/*<! number of simulated clusters */
	unsigned int *sim_cluster;	/*<! simulated cluster assignments */
	char const *sim_info_file;	/*<! info if data are simulated */
}; /* options */

int make_opt(options **opt);
void free_opt(options *opt);
int parse_opt(options *opt, int argc, const char **argv);
int process_arg_p(int argc, char const **argv, int *i, int j, options *opt);
const char *khaplotype_algorithm(int algorithm);

#endif /* HAPLOTYPE_OPTION_H */
