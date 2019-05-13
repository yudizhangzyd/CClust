/**
 * @file fqmorph.h
 * @author Karin S. Dorman
 *
 * Header file for fqmorph structs, extern functions and defines.
 */

#ifndef __H_FQMORPH__
#define __H_FQMORPH__

#include <stdio.h>
#include <stddef.h>

#include "fastq.h"
#include "error.h"

typedef struct _options options;
typedef struct _data data;

enum {
	FQMORPH_NO_ERROR,
	FQMORPH_AMBIGUOUS_NUCLEOTIDE	/* example only */
};

/**
 * Output format types.
 */
enum {
	FASTQ_FORMAT,
	FASTA_FORMAT,
	TEXT_FORMAT,	/* format used by k-modes */
	NUMBER_FORMATS
};

/**
 * Run options.
 */
struct _options {
	/* data */
	char const *fastq_file;	/*<! name of fastq input file */
	char const *partition_file;	/*<! partition file */

	/* running */
	int remove_ambiguous;	/*<! remove reads with ambiguous characters */
	int reverse_complement;	/*<! request reverse complement */
	int split_file;		/*<! use partition file to split */
	int output_clusters;	/*<! output cluster assignments */
	int append;		/*<! append output to existing file */
	int output_format;	/*<! format of output */
	int expected_errors;	/*<! expected number of errors */
	unsigned long seed;	/*<! random number seed [srand()] */
	size_t cut_start;	/*<! start of read region to cut */
	size_t cut_end;		/*<! end of read region to cut */
	size_t subsample;	/*<! randomly subsample this many reads */
	unsigned int nni_k;	/*<! k of k-nearest neighbor */

	char const *outfile;	/*<! output file */
}; /* options */

/**
 * Data.
 */
struct _data {
	unsigned int *cluster_id;	/*<! store contents of partition file */
	/* could add pointers to reads here */
}; /* data */

const char * fqmorph_error_message(int err_no);


#endif
