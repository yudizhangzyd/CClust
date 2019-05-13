#ifndef __IO_KMODES_H__
#define __IO_KMODES_H__

#include <stddef.h>	/* FILE, fprintf() */
#include <curses.h>	/* WINDOW, wprintw() */

#include "constants.h"

/* print vectors */
void fprint_data_ts(FILE *fp, data_t *v, size_t len, int width, int newline);
void wprint_data_ts(WINDOW *fp, data_t *v, size_t len, int width, int newline);

/* print matrices */
void fprint_data_t_matrix(FILE *fp, data_t **v, size_t n, size_t p, int width, int newline);
void wprint_data_t_matrix(WINDOW *fp, data_t **v, size_t n, size_t p, int width, int newline);

/* read vectors */
int fscan_data_ts(FILE *fp, data_t *v, size_t len);

#endif
