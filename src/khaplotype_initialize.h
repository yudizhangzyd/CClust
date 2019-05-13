#ifndef __KHAPLOTYPE_INITIALIZE_H__
#define __KHAPLOTYPE_INITIALIZE_H__

#include "khaplotype_data.h"
#include "khaplotype_model.h"
#include "khaplotype_option.h"
#include "init.h"

/**
 * Initialization methods are all those of k-modes plus one.
 */
enum {
	KHAPLOTYPE_INIT_FILTER_DATA = KMODES_INIT_NUMBER_METHODS + 1,
	KHAPLOTYPE_INIT_NUMBER_METHODS
};


int setup_initializer(data *dat, options *opt);
int initialize(data *dat, auxvar *aux_dat, options *opt, model *mod);
char const *khaplotype_init_method(int init_method);


#endif
