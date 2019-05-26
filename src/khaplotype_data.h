/**
 * @file khaplotype_data.h
 * @author Yudi Zhang
 */

#ifndef KHAPLOTYPE_H
#define KHAPLOTYPE_H

#include "fastq.h"
#include "constants.h"
#include "khaplotype_option.h"
#include "io_kmodes.h"
#include "io.h"
#include "hash.h"

typedef struct _data data;
typedef struct _auxvar auxvar;
typedef struct _outres outres;


/**
 * Data structure for clustering fastq data.
 */
struct _data {
	
	/* fastq data */
	fastq_data *fdata;
	data_t **qmat;			/*<! quality score matrix*/
	unsigned char n_quality;	/*<! number quality scores [min, max] */
	unsigned int max_read_length;	/*<! maximum observed read length */
	unsigned int min_read_length;	/*<! minimum observed read length */
	unsigned int *lengths;		/*<! length of reads */
	
	/* index */
	size_t *read_idx;		/*<! index array of all reads* */
	
	/* truth when data simulated */
	unsigned int *true_cluster_id;	/*<! true cluster assignments */
	unsigned int *true_cluster_size;/*<! true cluster sizes */
	
	data_t *data;			/*<! data */
	data_t **dmat;			/*<! data as matrix */
	unsigned int n_observations;	/*<! number of observations */
	unsigned int n_coordinates;	/*<! no. of coordinates (all categorical) */
	data_t *n_categories;		/*<! no. of categories per coordinate */
	data_t max_n_categories;	/*<! max. number categories per coordinate */
	data_t tot_n_categories;	/*<! total number of categories */
	data_t *categories;		/*<! total no. of the category */
	
	/* hash table */
	hash *seq_count; /*<! frequency of unique sequences (hash table) */
	unsigned int hash_length;  /*<! num of unique sequences in hash table */
	unsigned int hash_leng_exk;  /*<! num of unique sequences with abundance exceed k in hash table */
	unsigned int *seq_abundance; /*<! abundance exceed k in hash table */
	size_t *hash_idx_exk; /*<! index of unique sequences with abundance exceed k in hash table */
	size_t **hash_idx_exk_arr;
	double *llk_seq; /*<! sum of the log likelihood for unique sequences with abundance exceed k in hash table */
	unsigned int seed_count;
	size_t *min_index;
	
	/* mean expected error */
	double *exp_no_errs;		/*<! exp. no errors per read */
	double *mean_exp_no_errs;	/*<! mean exp. no errors per sequence */
	unsigned int *reads_uniq_id;	/*<! index of each read in unique sequence table */
	unsigned int *uniq_seq_count;	/*<! abundance of unique sequences */
	size_t *uniq_seq_idx;		/*<! unique sequences index in dmat */
	size_t **idx_arr;		/*<! index array under unique sequence */
	size_t *sel_uniq_seq_idx;
	size_t **sel_idx_arr;
	unsigned int *sel_uniq_seq_count;
	int n_selected_sequences;
	
	/* initialization */
	data_t **seeds;			/*<! chosen seeds */
	data_t **ini_seeds;		/*<! trial seeds */
	unsigned int *seed_idx;		/*<! seed indices */
	unsigned int *ini_seed_idx;	/*<! trial seed indices */
	unsigned int n_init;		/*<! number of initializations done */
	uint8_t use_ini;		/*<! using ini_* versions or not */
	
	/* current solution */
	double total;			/*<! current criterion */
	unsigned int *cluster_id;	/*<! cluster assignments */
	unsigned int *obsn_idx;		/*<! used if shuffling */
	double *criterion;		/*<! criterion */
	unsigned int *cluster_size;	/*<! cluster sizes */
	unsigned int iter;		/*<! iterations */
	
	/* best solution */
	double best_total;		/*<! total criterion */
	double best_rand;		/*<! if simulated */
	unsigned int *best_seed_idx;	/*<! seeding */
	data_t **best_modes;		/*<! estimated modes */
	double *best_criterion;		/*<! criterion */
	unsigned int *best_cluster_id;	/*<! cluster assignments */
	unsigned int *best_cluster_size;/*<! cluster sizes */
	unsigned int *best_obsn_idx;	/*<! used if shuffling */
	
	/* summary statistics */
	double seconds;			/*<! seconds used */
	double avg_cost, sd_cost;	/*<! minimum criterion */
	double avg_iter, sd_iter;	/*<! iterations to convergence */
	double avg_ar, sd_ar;		/*<! adjusted rand (truth known) */
	double avg_mi, sd_mi;		/*<! mutual information (truth known) */
	double avg_vi, sd_vi;		/*<! variance of info (truth known) */
	double avg_time, sd_time;	/*<! inits to target */
	unsigned int ntimes;		/*<! times hit target */
	
	/* internal use */
	double uncounted_seconds;
	double first_cost;
	double worst_cost;
	unsigned int max_iter;
	unsigned int ctime;
}; /* data */

struct _auxvar {
	/* auxiliary variables */
	double **prob_t;		/*<! probability of no error: log(1-p_ij) */
	double **prob_f;		/*<! probability of error: log(p_ij/3) */
	double ***e_kjN;		/*<! likelihood: nuc. N at site j in clus. k */
	double **v_ik;			/*<! likelihood: obs. i in cluster k */
	unsigned int *last_assign;	/*<! previous assignment of observations */
	double *last_cost;		/*<! previous cost */
	double **last_vik;		/*<! previous vik */
	double ***last_ekjn;		/*<! previous ekjn */
	unsigned int **last_cent;	/*<! previous centers */
	unsigned int *all_index;	/*<! when all clusters live */
	unsigned int *some_index;	/*<! ordered indices of live clusters */
	unsigned int *index;		/*<! pointer to current live cluster index */
	unsigned int *index_in_cluster; /*<! pointer to store the indices of obervation */
	unsigned int *live;			/*<! indicate if it is a live set */
	double *llk_diff;		/*<! difference of the likelihood in each cluster*/
	unsigned int **idx_in_cluster;	/*<! index of each observation is each cluster */
	unsigned int *idx_position;	/*<! position of each observation is cluster */
	unsigned int **changed_center;
	unsigned int *cent_j_len;
};

struct _outres {
	unsigned int n_observations;	/*<! number of observations */
	unsigned int n_coordinates;	/*<! no. of coordinates (all categorical) */
	/* best solution */
	double best_total;		/*<! total criterion */
	unsigned int *best_modes;		/*<! estimated modes */
	double *best_criterion;		/*<! criterion */
	unsigned int *best_cluster_id;	/*<! cluster assignments */
	unsigned int *best_cluster_size;/*<! cluster sizes */
	unsigned int *data;
};

int make_hap_data(data **data);
int sync_state_khaplotype_data(data *dat, options *opt);
int read_hap_data(data *dat, options *opt);
int finish_make_hap_data(data *dat, options *opt);
void free_hap_data(data *dat);
void free_aux_var(auxvar *aux);
int make_aux_var(auxvar **aux);
int sync_state_aux_var(data *dat, auxvar *aux, options *opt);
int write_solution(data *dat, options *opt, FILE **in_fps);
void sample_better(unsigned int N, unsigned int n, unsigned int *idx);
int read_fsa(const char *filename, int n, int p, data_t **true_seed);
int make_seeds(data *dat, options *opt);
int shuffle_data_hap(data *dat, options *opt);
int build_hash(data *dat, options *opt);
int filter_haplotypes(data *dat, options *opt);
void compute_pij(void *d, double **prob_t, double **prob_f);
int make_res(outres **out, data *dat, options *opt);
void free_res(outres *out);

#endif /* HAPLOTYPE_H */
