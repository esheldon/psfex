/*
Some wrappers for fitsio that automatically report errors
and simplify some tasks

Only routines where error can only be determined by status require
status keywords.

Also provide functions for most cases, e.g. reading a key
*/
#include <fitsio.h>

#define FT_DEFVAL -9999

fitsfile *ft_open(const char *filename);
fitsfile *ft_close(fitsfile *fits);
double ft_read_key_dbl(fitsfile *fits, const char *name, int *status);
long ft_read_key_lng(fitsfile *fits, const char *name, int *status);
int ft_get_colnum(fitsfile *fits, const char *colname, *status);
long ft_get_array_col_size(fitsfile *fits, const char *colname, int *status);

char *ft_read_str_row_byname(fitsfile *fits,
                             const char* colname,
                             long row,
                             int *status);

static double fits_read_dbl_row_byname(fitsfile *fits,
                                       const char* colname, 
                                       long row,
                                       int *status);
