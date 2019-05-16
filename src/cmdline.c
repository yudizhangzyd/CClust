/**
 * @file cmdline.c
 * @author Karin S. Dorman, kdorman@iastate.edu
 * @date Tue Mar 20 18:10:33 CDT 2012
 *
 * Some functions for processing a command line.  Note, the programmer must
 * provide the definition of a function void fprintf_usage(FILE *, const char
 * *, void *) in order to use this code.  This function, as its name suggests,
 * prints help on command-line usage for the user.
 *
 * Note: The requirement used to be that the user defined print_usage(const
 * char *).  To compile older code under that obsolete system, use command-line
 * option -DOLD_USAGE with gcc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "cmdline.h"

/**
 * Print error message, usage information, and return a failure code.
 *
 * @param argv list of comman-line arguments as (unmutable) string constants
 * @param i index of argument that could not be parsed
 * @param obj void pointer to additional usage information
 * @return (simple) error code
 */
int usage_error(const char **argv, int i, void *obj)
{
	fprintf(stderr, "ERROR -- incorrect command line usage\n");
	if (i >= 0)
		fprintf(stderr, "Could not parse command line option '%s'.\n",
			argv[i]);
	fprintf(stderr, "See correct usage below.\n");
#ifdef OLD_USAGE	/* define to compile with older source code */
	print_usage(argv[0]);
#else
	fprint_usage(stderr, argv[0], obj);
#endif
	fprintf(stderr, "\nERROR -- incorrect command line usage\n");
	if (i >= 1)
		fprintf(stderr, "Could not parse command line option '%s'.\n",
			argv[i]);
	return EXIT_FAILURE;
} /* usage_error */

/**
 * Read an integer from the next argument in the command-line argument
 * list.  This code avoids some repetition in main()
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed integer or errno is error
 */
int read_int(int argc, const char **argv, int i, void *obj)
{
	int n;
	char *ret_ptr;

	/* check for end of argument list; assignment to errno BEFORE return */
	if (i == argc)
		return (errno = usage_error(argv, i-1, obj));
	/* It turns out errno is not set to EINVAL (in C99) when the read fails
	 * as I said in class.  The only way to detect failure to read a
	 * numeric value is to use the second argument.   The second argument
	 * is set to point to the char immediately following the read number.
	 * If it is the same as the address in the first pointer passed, then
	 * you know the read was unsuccessful.  Note the need to pass a pointer
	 * to a pointer in order to update the value (address) stored in the
	 * pointer variable. */
	/* n = strtol(argv[i], NULL, 0); */
	errno = 0;
	n = strtol(argv[i], &ret_ptr, 0);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));
	return n;
} /* read_int */

/**
 * Read an positive integer from the next argument in the command-line argument
 * list.  This code avoids some repetition in main()
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed unsigned integer or errno is error
 */
unsigned int read_uint(int argc, const char **argv, int i, void *obj)
{
	unsigned int n;
	char *ret_ptr;

	/* check for end of argument list; assignment to errno BEFORE return */
	if (i == argc)
		return (errno = usage_error(argv, i-1, obj));
	/* It turns out errno is not set to EINVAL (in C99) when the read fails
	 * as I said in class.  The only way to detect failure to read a
	 * numeric value is to use the second argument.   The second argument
	 * is set to point to the char immediately following the read number.
	 * If it is the same as the address in the first pointer passed, then
	 * you know the read was unsuccessful.  Note the need to pass a pointer
	 * to a pointer in order to update the value (address) stored in the
	 * pointer variable. */
	/* n = strtol(argv[i], NULL, 0); */
	errno = 0;
	n = strtoul(argv[i], &ret_ptr, 0);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));
	return n;
} /* read_uint */

/**
 * Read a long from the next argument on the command-line.
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed long or errno is error
 */
long read_long(int argc, const char **argv, int i, void *obj)
{
	long n;
	char *ret_ptr;

	if (i == argc)
		return (errno = usage_error(argv, i-1, obj));

	errno = 0;
	n = strtol(argv[i], &ret_ptr, 0);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));

	return n;
} /* read_long */

/**
 * Read an unsigned long from the next argument on the command-line.
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed long or errno is error
 */
unsigned long read_ulong(int argc, const char **argv, int i, void *obj)
{
	unsigned long n;
	char *ret_ptr;

	if (i == argc)
		return (errno = usage_error(argv, i-1, obj));

	errno = 0;
	n = strtoul(argv[i], &ret_ptr, 0);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));

	return n;
} /* read_ulong */

/**
 * Read a double from the next argument in the command-line argument
 * list.  This code avoids repetition in main().
 *
 * @param argc the number of arguments
 * @param argv the list of command-line arguments strings
 * @param i the argument we are to parse at int
 * @param obj void pointer to additional usage information
 * @return the parsed integer or errno is error
 */
double read_cmdline_double(int argc, const char **argv, int i, void *obj)
{
	double d;
	char *ret_ptr;

	if (i==argc)
		return (errno = usage_error(argv, i-1, obj));

	errno = 0;
	d = strtod(argv[i], &ret_ptr);
	if (errno || ret_ptr == argv[i])
		return (errno = usage_error(argv, i, obj));

	return d;
} /* read_cmdline_double */

/**
 * Read a char from the next argument on the command-line.
 *
 * @param argc number of arguments on command-line
 * @param argv list of command-line argument strings
 * @param i the argument we are to parse as char
 * @param obj void pointer to additional usage information
 * @return the parsed char or errno is error
 */
char read_char(int argc, const char **argv, int i, void *obj)
{
	if (i==argc)
		return (errno = usage_error(argv, i-1, obj));

	return argv[i][0];
} /* read_char */

/**
 * Read a list of doubles from the next arguments.
 */
unsigned int read_cmdline_doubles(unsigned int argc, const char **argv,
	unsigned int i, double **dret, void *obj)
{
	unsigned int n = 0;
	while (i + n < argc) {
		if (argv[i + n][0] == '-')
			break;
		else
			n++;
	}
	if (!n)
		return (errno = usage_error(argv, i-1, obj));
	if (!dret)
		return n;

	*dret = malloc(n * sizeof **dret);
	if (*dret)
		for (unsigned int j = 0; j < n; ++j)
			(*dret)[j] = read_cmdline_double(argc, argv, i++, obj);
	return n;
} /* read_cmdline_doubles */

/**
 * Read a list of unsigned long from the next arguments.
 */
unsigned int read_ulongs(unsigned int argc, const char **argv, unsigned int i,
	unsigned long **sret, void *obj)
{
	unsigned int n = 0;
	while (i + n < argc) {
		if (argv[i + n][0] == '-')
			break;
		else
			n++;
	}
	if (!n)
		return (errno = usage_error(argv, i-1, obj));
	if (!sret)
		return n;

	*sret = malloc(n * sizeof **sret);
	if (*sret)
		for (unsigned int j = 0; j < n; ++j)
			(*sret)[j] = read_ulong(argc, argv, i++, obj);

	return n;
} /* read_ulongs */

/**
 * Read a list of unsigned ints from the next arguments.
 */
unsigned int read_cmdline_uints(unsigned int argc, const char **argv,
	unsigned int i, unsigned int **uiret, void *obj)
{
	unsigned int n = 0;
	while (i + n < argc) {
		if (argv[i + n][0] == '-')
			break;
		else
			n++;
	}
	if (!n)
		return (errno = usage_error(argv, i-1, obj));
	if (!uiret)
		return n;

	*uiret = malloc(n * sizeof **uiret);
	if (*uiret)
		for (unsigned int j = 0; j < n; ++j)
			(*uiret)[j] = read_uint(argc, argv, i++, obj);

	return n;
} /* read_cmdline_uints */

/**
 * Read a list of strings from the next arguments.
 */
unsigned int read_cmdline_strings(unsigned int argc, const char **argv,
	unsigned int i, char const ***sret, void *obj)
{
	unsigned int n = 0;
	while (n + i < argc) {
		if (argv[i + n][0] == '-')
			break;
		else
			++n;
	}
	if (!n)
		return (errno = usage_error(argv, i-1, obj));

	if (!sret)
		return n;

	*sret = malloc(n * sizeof **sret);
	if (*sret)
		for (unsigned int j = 0; j < n; ++j)
			(*sret)[j] = argv[i++];

	return n;
} /* read_cmdline_strings */

int is_numeric(const char *argv)
{
	size_t len = strlen(argv);
	for (size_t i = 0; i < len; ++i)
		if ((argv[i] < 48 || argv[i] > 57) && argv[i] != '.'
			&& argv[i] != '-' && argv[i] != 'e')
			return 0;
	return 1;
}/* is_numeric */


/**
 * print command line usage
 *
 * @param fp file
 * @param cmdname command line
 * @param obj option struct
 */
void fprint_usage(FILE *fp, const char *cmdname, void *obj)
{
	size_t start = strlen(cmdname) - 1;
	
	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;
	
	for (size_t i = start; i < strlen(cmdname); ++i) fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - cluster observations with categorical predictors\n", &cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s [-r <rulong> -n <nuint> -a|-e|-l|-w|-d -i <istr>|<iuint1...k> -o <ofile>] -k <kuint> -f <ffile> ...\n", &cmdname[start]);
	fprintf(fp, "\nDESCRIPTION\n\t%s clusters observations found in file <ffile> into <kuint> clusters.  It randomly initialize <iuint> times after setting random number seed <sulong>.\n", &cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t-k <kuint>\n\t\tSet the desired number of clusters K.\n");
	fprintf(fp, "\t-x <kuint>\n\t\tSet the cut-off value for mean expected error(default is 2.0).\n");
	fprintf(fp, "\t-z <kuint>\n\t\tSet the cut-off value for abundance in the hash table(default is 3).\n");
	fprintf(fp, "\t-y <kuint>\n\t\tSet the cut-off value for the proportion of BIC(default is 0.005).\n");
	fprintf(fp, "\t-b\n\t\tTest on categorical data without quality score [NOTICE: Hartigan and Wong is slow compared to other k-modes method!].\n");
	fprintf(fp, "\t-l\n\t\tRun Lloyd's algorithm.\n");
	fprintf(fp, "\t-e\n\t\tRun Lloyd's efficient algotithm for haplotype.\n");
	fprintf(fp, "\t-a\n\t\tRun Macqueen(Huang)'s algorithm.\n");
	fprintf(fp, "\t-w\n\t\tRun Hartigan and Wong algorithm.\n");
	fprintf(fp, "\t-d\n\t\tRun Efficient Hartigan and Wong algorithm.\n");
	fprintf(fp, "\t--shuffle\n\t\tShuffle the data input order on each initialization.\n");
	fprintf(fp, "\t-s <sn> <sp> <sc> <st1> <st2>\n\t\tSimulation size (observations by coordinates by categories) and times.\n");
	fprintf(fp, "\t-f <ffile> \n\t\tSet the input filename.\n");
	fprintf(fp, "\t-o <ofile> \n\t\tSet the output filename.\n");
	fprintf(fp, "\t-i <istr>\n\t\tSet initialization method (rnd|filter).\n");
	fprintf(fp, "\t\t'rnd' randomly selects seeds.\n");
	fprintf(fp, "\t\t'filter' filtering by hash table and expected mean error to select seeds.\n");
	fprintf(fp, "\t-r <rulong>\n\t\tSet random number seed. [OPTIONAL]\n");
	fprintf(fp, "\t-h\n\t\tThis help.\n");
	fprintf(fp, "\n");
	for (size_t i = start; i < strlen(cmdname); ++i) fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
} /* fprint_usage */
