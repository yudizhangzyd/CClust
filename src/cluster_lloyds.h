/**
 * @file cluster_lloyds.h
 * @author Yudi Zhang
 *
 * Header file for culster_lloyds
 *
 */

#ifndef CLUSTER_LLOYDS_H
#define CLUSTER_LLOYDS_H

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

#define LOG3 1.098612288668109782108

typedef int (*lloyd_step1)(void *, void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *, unsigned int *, int);
typedef void (*lloyd_step2)(void *, void *, data_t **, unsigned int, unsigned int, unsigned int, unsigned int *);
typedef double (*compute_rule)(void *, void *, data_t **, unsigned int *ic, double *, unsigned int K, unsigned int n, unsigned int p);

data_t *find_unique(data_t *mat, unsigned int n, data_t *n_unique);
int compare(const void *a, const void *b);

double compute_criterion_hap(void *d, void *auxd, data_t **seeds, unsigned int *ic, double *criterion, unsigned int K, unsigned int n, unsigned int p);
int fastq_lloyds_step1(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *nclass, unsigned int *ic1, int Is_Init);
void fastq_lloyds_step2(void *d, void *auxd, data_t **seeds, unsigned int K, unsigned int n, unsigned int p, unsigned int *ic1);

double cluster_lloyds(void *d, void *auxd, data_t **seeds, unsigned int *nclass, unsigned int *ic1, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *iter, lloyd_step1 step1, lloyd_step2 step2, compute_rule compute_sum_cost);

#endif /* cluster_lloyds_h */
