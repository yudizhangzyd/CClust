/*
 *  @file init.h
 */

#ifndef __KMODES_INIT_H__
#define __KMODES_INIT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>

#include "math.h"
#include "constants.h"
#include "error.h"
#include "sample.h"
#include "kmodes.h"

/**
 * K-modes initialization methods.
 */
enum {
	KMODES_INIT_RANDOM_SEEDS,	/*!< klaR() */
	KMODES_INIT_RANDOM_FROM_PARTITION,	/*!< random sampling from true clusters */
	KMODES_INIT_RANDOM_FROM_SET,	/*!< random sampling from set */
	KMODES_INIT_NUMBER_RANDOM_METHODS
};


const char *kcluster_init_method(int init_method);

int compare_data(data_t **x, unsigned int i, unsigned int j, unsigned int n);
int kcluster_init(data_t **x, unsigned int n, unsigned int p, unsigned int k, unsigned int k1, data_t **seeds, unsigned int *seed_idx, int method);
int kmodes_init_random_from_partition(data_t **x, unsigned int n, unsigned int p, unsigned int K, data_t **seeds, unsigned int *sidx, unsigned int *id);
int kmodes_init_random_from_set(unsigned int K, unsigned int p, unsigned int n_seedset, data_t **seeds, data_t **seedset);

#endif /* kmodes_init.h */
