/**
 * @file kmodes.h
 *
 * [TODO] In the process of separating initializaton methods from k-modes code 
 * (see kmodes_initialization.[ch]), but it will require getting rid of the 
 * global variables that k-modes code uses.  The better solution is to use
 * Yudi's code for k-modes once it is confirmed to be faster.
 *
 * WARNING: user must define data type data_t before including this file.
 */
#ifndef KMODES_H
#define KMODES_H

#include <stddef.h>

#include "constants.h"
#include "init.h"

/**
 * Default maximum number of K-modes iterations.
 * The user can specify any number of maximum iterations when calling
 * #kmodes(), but for consistency across the code, this define can be
 * used.
 */
#define MAX_KMODES_ITERATIONS	1000

/**
 * Errors produced by kmodes().
 */
enum {
	KMODES_NO_ERROR,		/*!< error-free condition */
	KMODES_EXCEED_ITER_WARNING,	/*!< K-modes exceeded max iterations */
	KMODES_NO_WARNINGS,		/*!< no. non-lethal warnings */
	KMODES_NULL_CLUSTER_ERROR,	/*!< K-modes produced null cluster */
	KMODES_CALLER_INPUT_ERROR,	/*!< invalid caller input */
	KMODES_MEMORY_ERROR,		/*!< memory allocation error */
	KMODES_INVALID_INITIALIZATION_METHOD,
	KMODES_INTERNAL_ERROR,
	KMODES_NUMBER_ERRORS		/*!< total number of kmodes error codes */
};

/**
 * K-modes algorithms.
 */
enum {
	KMODES_HUANG,			/*!< Huang 1997 */
	KMODES_LLOYD,			/*!< Lloyd 1982 */
	KMODES_HARTIGAN_WONG,		/*!< Hartigan & Wong 1979 */
	KMODES_NUMBER_ALGORITHMS
};

/**
 * K-modes variations.
 */
enum {
	KMODES_UPDATE_MODE,		/*!< H97 & HW79 update modes... */
	KMODES_NUMBER_VARIATIONS	/*   ...during initialization */
};

/**
 * K-modes initial states.  Each initialization method must produce
 * an initial state of seeds or a partition.
 */
enum {
	KMODES_INIT_STATE_SEEDS,	/*!< seeds */
	KMODES_INIT_STATE_PARTITION,	/*!< partition */
	KMODES_INIT_STATE_NUMBER
};

/**
 * K-modes weighting methods
 */
enum {
	KMODES_NO_WEIGHTING,
	KMODES_NUMBER_WEIGHTING_METHODS
};

/**
 * format of the pointer of the index
 */
enum {
	UNSINGED_INT,
	_SIZE_T_P
};

/**
 * K-modes initialization methods.
 */
enum {
	KMODES_INIT_RANDOM_SEEDS_AUX,	/*!< klaR() */
	KMODES_INIT_RANDOM_FROM_PARTITION_AUX,	/*!< random sampling from true clusters */
	KMODES_INIT_RANDOM_FROM_SET_AUX,	/*!< random sampling from set */
	KMODES_INIT_NUMBER_RANDOM_METHODS_AUX,
	KMODES_INIT_H97_RANDOM,		/*!< Huang1997, randomized version */
	KMODES_INIT_HD17,		/*!< Huang1997 interpretted by Python author de Vos */
	KMODES_INIT_CLB09_RANDOM,	/*!< Cao2009, randomized version */
	KMODES_INIT_AV07,		/*!< K-means++ adapted */
	KMODES_INIT_AV07_GREEDY,	/*!< K-means++ greedy adapted */
	KMODES_INIT_USER_SEEDS,		/*!< user provides K seeds */
	KMODES_INIT_H97,		/*!< Huang1997 */
	KMODES_INIT_CLB09,		/*!< Cao2009 */
	KMODES_INIT_NUMBER_METHODS
};

/**
 * Small object with k-modes options for run variations.  Many more options are
 * found in run_kmodes.c:options.
 */
typedef struct _kmodes_options kmodes_options;
struct _kmodes_options {
	int weighted;		/*!< use Huang's weighted version */
	int init_update;	/*!< update modes during initialization */
	int use_qtran;		/*!< use quick transfer in Hartigan & Wong */
	int use_hartigan;	/*!< use Hartigan's transfer decision rule */
};

int make_kmodes_options(kmodes_options **opt);

size_t allocate_and_compute_category_counts(data_t **x, unsigned int n, unsigned int p, unsigned int K, int wgt);
int set_nkjc(unsigned int K, unsigned int p, size_t l);
void compute_nkjc(size_t ***nkjc, data_t **x, unsigned int *ic, unsigned int n, unsigned int p);

/* k-modes algorithm */
double kmodes_hw(data_t **a, unsigned int m, unsigned int n, data_t **c, unsigned int k, unsigned int *ic1, unsigned int *nc, unsigned int miter, double *wss, int *ifault, unsigned int *iter, kmodes_options *opt);//int weight, int init_update, int use_qtran);
double kmodes_huang(data_t **x, unsigned int n, unsigned int p, data_t **seeds, unsigned int K, unsigned int *ic1, unsigned int *nclass, unsigned int miter, double *wss, int *ifault, unsigned int *iter, kmodes_options *opt);//int weight, int init_update);
double kmodes_lloyd(data_t **x, data_t **seeds, unsigned int *nclass, unsigned int *ic1, unsigned int n, unsigned int p, unsigned int K, unsigned int miter, double *wss, int *ifault, unsigned int *iter, int weight);


const char *kmodes_algorithm(int algorithm);

const char *kmodes_init_method(int init_method);
int kmodes_init_from_partition(data_t **x, unsigned int n, unsigned int p, unsigned int k, int weight, data_t **seeds, unsigned int *ic);
int kmodes_init (data_t **x, unsigned int n, unsigned int p, unsigned int k, unsigned int k1, data_t **seeds, unsigned int *sidx, int method, int weight);

/*
void optra(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, unsigned int *ic1,
        unsigned int *ic2, SIZE_T *nc, SIZE_T *live, SIZE_T *cd, double *d, SIZE_T *indx,
	SIZE_T *nj, SIZE_T *nt, SIZE_T *nkt);
void qtran(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, unsigned int *ic1,
        unsigned int *ic2, SIZE_T *nc, SIZE_T *live, double *d, SIZE_T *indx,
	SIZE_T *nj, SIZE_T *nt, SIZE_T *nkt);
*/

const char *kmodes_error(int err);
void free_kmodes();	/* use when rerunning with altered data */
void reset_k();		/* use when rerunning with altered K */

#endif /* KMODES_H */
