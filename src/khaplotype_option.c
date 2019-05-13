/**
 * @file khaplotype_option.c
 * @author Yudi Zhang
 */

#include <unistd.h>

#include "khaplotype.h"
#include "khaplotype_option.h"
#include "khaplotype_initialize.h"
#include "kmodes.h"
#include "error.h"
#include "math.h"
#include "fastq.h"
#include "io.h"
#include "io_kmodes.h"
#include "cmdline.h"

const char *khaplotype_algorithm(int algorithm)
{
	if (algorithm == FASTQ_MACQUEEN)
		return "MacQueen 1997";
	else if (algorithm == FASTQ_LLOYDS)
		return "Lloyd 1982";
	else if (algorithm == FASTQ_LLOYDS_EFFICIENT)
		return "Efficient Lloyd 1982";
	else if (algorithm == FASTQ_HW)
		return "Hartigan Wong 1979";
	else if (algorithm == FASTQ_HW_EFFICIENT)
		return "Efficient Hartigan Wong 1979";
	else
		return NULL;
} /* khaplotype_algorithm */

int make_opt(options **opt)
{
	*opt = malloc(sizeof **opt);
	
	if (*opt == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "options object");
	
	(*opt)->K = 0;    /* invalid value */
	(*opt)->true_column = UINT_MAX;
	(*opt)->run_with_quals = 1;
	(*opt)->true_cluster = NULL;
	(*opt)->true_K = 0;
	(*opt)->hash_k = 3;
	(*opt)->cut_k = 0.005;
  	(*opt)->cutoff_mean_exp_err = 0.8;
	(*opt)->filter_haplotypes = 1;
	(*opt)->run_with_hash = 1;
	(*opt)->n_init = 1; /* Used to test the running time */
	(*opt)->n_inner_init = 1;
	(*opt)->n_max_iter = 10000;
	(*opt)->info = QUIET;
	(*opt)->subtract_one = 0;
	(*opt)->shuffle = 0;
	(*opt)->khaplotype_algorithm = FASTQ_LLOYDS_EFFICIENT;
	(*opt)->init_method = KMODES_INIT_RANDOM_SEEDS;
	(*opt)->seed = 0;
	(*opt)->seed_idx = NULL;
	(*opt)->n_sd_idx = 0;
	(*opt)->n_seed_set = 0;
	(*opt)->datafile = NULL;
	(*opt)->soln_file = NULL;
	(*opt)->quiet = MINIMAL;
	(*opt)->sim_K = 0;
	(*opt)->sim_cluster = NULL;
	(*opt)->sim_info_file = NULL;
	
	return NO_ERROR;
} /* make_options */

void free_opt(options *opt)
{
	if (opt) {
		if (opt->true_cluster) {
			free(opt->true_cluster);
			opt->true_cluster = NULL;
		} else if (opt->sim_cluster) {
			free(opt->sim_cluster);
			opt->sim_cluster = NULL;
		}
		/* opt->seed_idx free'd via dat->seed_idx */
		free(opt);
	}
} /* free_options */

int parse_opt(options *opt, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;
	
	for (i = 1; i < argc; i++) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];
		switch(a) {
			case 'b':
				opt->run_with_quals = 0;
				break;
			case 'v':
				opt->filter_haplotypes = 0;
				break;
			case 'k':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
					opt->K = read_uint(argc, argv, ++i,
							   (void *)opt);
					if (errno)
						goto CMDLINE_ERROR;
					debug_msg(MINIMAL, opt->quiet,
						  "K = %u.\n", opt->K);
				break;
			case 'l':
				opt->khaplotype_algorithm = FASTQ_LLOYDS;
				debug_msg(QUIET, opt->quiet,
					  "Using Lloyd's algorithm.\n");
				break;
			case 'e':
				opt->khaplotype_algorithm = FASTQ_LLOYDS_EFFICIENT;
				debug_msg(QUIET, opt->quiet, "Using efficient "
					"Lloyd's algorithm.\n");
				break;
			case 'w':
				opt->khaplotype_algorithm = FASTQ_HW;
				debug_msg(QUIET, opt->quiet,
					  "Using Hartigan and Wong algorithm.\n");
				break;
			case 'd':
				opt->khaplotype_algorithm = FASTQ_HW_EFFICIENT;
				debug_msg(QUIET, opt->quiet, "Using Hartigan "
					"and Wong algorithm (efficient).\n");
				break;
			case 'a':
				opt->khaplotype_algorithm = FASTQ_MACQUEEN;
				debug_msg(QUIET, opt->quiet,
					  "Using Macqueen's algorithm.\n");
				break;
			case 'f':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->datafile = argv[++i];
				debug_msg(QUIET, opt->quiet, "Data file = %s\n",
					  opt->datafile);
				break;
			case 'o':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->soln_file = argv[++i];
				break;
			case 'i':
				if (i + 1 == argc || argv[i + 1][0] == '-')
					goto CMDLINE_ERROR;
				if (!strcmp(argv[i+1], "rnd"))
					opt->init_method = KMODES_INIT_RANDOM_SEEDS;
				else if (!strcmp(argv[i+1], "filter"))
					opt->init_method = KHAPLOTYPE_INIT_FILTER_DATA;
				debug_msg(QUIET, opt->quiet, "Using %s "
					  "initialization.\n",
					  khaplotype_init_method(opt->init_method));
				++i;
				break;
			case 'n':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->n_init = read_uint(argc, argv, ++i,
							(void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				debug_msg(MINIMAL, opt->quiet,
					  "Initializations = %u\n", opt->n_init);
				break;
			case 'y':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->cut_k = read_cmdline_double(argc, argv, ++i,
								 (void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				debug_msg(MINIMAL, opt->quiet,
					  "cut_K = %lf.\n", opt->cut_k);
				break;
			case 'z':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->hash_k = read_uint(argc, argv, ++i,
							(void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				debug_msg(MINIMAL, opt->quiet,
					  "hash_K = %u.\n", opt->hash_k);
				break;
			case 'x':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->cutoff_mean_exp_err =
				read_cmdline_double(argc, argv,
						    ++i, (void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				debug_msg(MINIMAL, opt->quiet,
					  "mean_expected_error = %lf.\n",
					  opt->cutoff_mean_exp_err);
				break;
			case 'r':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				
				opt->seed = read_uint(argc, argv, ++i,
						      (void *)opt);
				srand(opt->seed);
				debug_msg(MINIMAL, opt->quiet, "Seed: "
					  "%lu\n", opt->seed);
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 's':    /*-s <sn> <sp> <sc> <st1> <st2>*/
				/* no appropriate arguments */
				if (!strcmp(&argv[i][j], "shuffle")) {
					opt->shuffle = 1;
					debug_msg(MINIMAL, opt->quiet, "Data "
						  "will be shuffled.\n");
					break;
				}
				if (i + 1 < argc
					&& access(argv[i+1], F_OK) != -1) {
					opt->sim_info_file = argv[++i];
					break;
				}
				if (errno)
					goto CMDLINE_ERROR;
				
				break;
			case 'h':
				fprint_usage(stderr, argv[0], opt);
				free_opt(opt);
				exit(EXIT_SUCCESS);
				break;
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}
	
	if (opt->K == 1)
		opt->n_init = 1;
	
	if (!opt->shuffle && opt->n_init > 1)
		mmessage(WARNING_MSG, NO_ERROR, "Failing to shuffle the data "
			 "can produce input order-determined behavior.\n");
	
	return err;
	
CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */
