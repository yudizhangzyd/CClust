/**
 * @file khaplotype.h
 * @author Yudi Zhang
 *
 * Header file for khaplotype code.
 *
 */

#ifndef __KHAPLOTYPE_H__
#define __KHAPLOTYPE_H__

#include <float.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <curses.h>
#include <math.h>

#define USE_CURSES

#include "constants.h"
#include "cluster_lloyds.h"

/**
 * Errors produced by khaplotype().
 */
enum {
	KHAPLOTYPE_NO_ERROR,		/*!< error-free condition */
	KHAPLOTYPE_EXCEED_ITER_WARNING,	/*!< K-modes exceeded max iterations */
	KHAPLOTYPE_NO_WARNINGS,		/*!< no. non-lethal warnings */
	KHAPLOTYPE_NULL_CLUSTER_ERROR,	/*!< K-modes produced null cluster */
	KHAPLOTYPE_CALLER_INPUT_ERROR,	/*!< invalid caller input */
	KHAPLOTYPE_MEMORY_ERROR,	/*!< memory allocation error */
	KHAPLOTYPE_INVALID_INITIALIZATION_METHOD,
	KHAPLOTYPE_INTERNAL_ERROR,
	KHAPLOTYPE_NUMBER_ERRORS	/*!< total number of kmodes error codes */
};

typedef void (*modified_lloyd_step2)(void *, void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *, int);
typedef void (*macq_ini)(void *, void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *, unsigned int *);
typedef int (*macq_iter)(void *, void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *, unsigned int *);
typedef double (*compute_crit)(void *, unsigned int *, double *, unsigned int, unsigned int);
typedef void (*hw_init)(void *,  void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *, unsigned int *, double *);
typedef void (*hw_iter)(void *, void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *, unsigned int *, double *, unsigned int, int, unsigned int *);
typedef void (*hw_fast_init)(void *, void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *, unsigned int *);
typedef void (*hw_fast_iter)(void *, void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *, unsigned int *, unsigned int, unsigned int *);

int fastq_lloyds_efficient_step1(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic1, int Init);
void fastq_lloyds_efficient_step2(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *ic1, int Init);

int fastq_macqueen_iter(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic1);
void fastq_macqueen_ini(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic1);

void hw_init_default(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic, double *cost);
void hw_fast_init_default(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic);

void hw_iter_default(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic, double *cost, unsigned int i, int Is_Live, unsigned int *indx);
void hw_fast_iter_default(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic, unsigned int i, unsigned int *indx);

void seed_llk(void *d, void *auxd, unsigned int length, unsigned int *abundance, size_t *unique_seq, size_t ** seq_arr);

double compute_cost(void *d, unsigned int *ic, double *criterion, unsigned int K, unsigned int n);

double cluster_lloyds_efficient(void *d, void *auxd, data_t **seeds, unsigned int *nclass, unsigned int *ic1, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *iter, lloyd_step1 step1, modified_lloyd_step2 step2, compute_crit compute_sum_cost);

double cluster_macqueen(void *d, void *auxd, unsigned int n, unsigned int p, data_t **seeds, unsigned int K, unsigned int *ic1, unsigned int *nclass, unsigned int max_iter, double *cost, int *ifault, unsigned int *iter, macq_ini ini, macq_iter itr, compute_crit compute_sum_cost);

double cluster_hw(void *d, void *auxd, data_t **seeds, unsigned int *nclass, unsigned int *ic, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *niter, hw_init init, hw_iter iter);
double cluster_hw_fast(void *d, void *auxd, data_t **seeds, unsigned int *nclass, unsigned int *ic, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *niter, hw_fast_init init, hw_fast_iter iter);

char const *khaplotype_error(int err);

#endif /* khaplotype.h */
