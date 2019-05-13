#ifndef KHAPLOTYPE_MODEL_H
#define KHAPLOTYPE_MODEL_H

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <curses.h>

#include "khaplotype_data.h"
#include "khaplotype_option.h"
#include "statistics.h"
#include "array.h"

typedef struct _model model;

struct _model {

	unsigned char *haplotypes;	/*<! haplotypes */
	double best_ll;			/*<! best log likelihood so far */
	
	double *samp_prob;
	unsigned int ind;
	
	unsigned int n_param;		/*<! number of parameters */
	
	/* choose a place to start randomization */
	double JC_ll;			/*<! log likelihood of haplotypes */
	double aic;
	double bic;
	double previous_bic;
	double percent_change;
	
	double *distance;		/*<! Kx1 distances from haplotypes to the ancestor */
	unsigned char * est_ancestor;	/*<! estimated ancestor */
};

int make_model(model **mod, data *dat);
int realloc_model(model *mod, data *dat);
void free_model(model *mod);
int modified_ic(unsigned char* hap, unsigned char *est_anc, double *distance, double best_ll, unsigned int K, double *JC_ll, double *n_aic, double *n_bic, unsigned int n_param, unsigned int max_read_length,size_t sample_size);
int m_JC69(unsigned char * hap, unsigned char * anc, double *dist, unsigned int K, unsigned int len);
double e_JC69(unsigned char * hap, unsigned char * anc, double *dist, unsigned int K, unsigned int len);
int sample_from_categorical_distribution(data *dat, model *mod);
int detect_true (data_t **optimal_seq, data_t **true_seq, unsigned int K, unsigned int p);
#endif /* khaplotype_model_h */
