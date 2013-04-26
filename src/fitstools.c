#include <stdio.h>
#include <stdlib.h>
#include "fitstools.h"

fitsfile *ft_open(const char *filename)
{
    int status=0;
    fitsfile *fits=NULL;
    if (fits_open_file(&fits, filename, READONLY, &status)) {
        fits_report_error(stderr,status);
        return NULL;
    }
    return fits;
}
fitsfile *ft_close(fitsfile *fits)
{
    int status=0;
    if (fits_close_file(fits, &status)) {
        fits_report_error(stderr,status);
    }
    return NULL;
}

double ft_read_key_dbl(fitsfile *fits, const char *name, int *status)
{
    double val=FT_DEFVAL;
    if (fits_read_key_dbl(fits,(char*)name,&val,NULL, status)) {
        fits_report_error(stderr,*status);
    }
    return val;
}
long ft_read_key_lng(fitsfile *fits, const char *name, int *status)
{
    long val=FT_DEFVAL;
    if (fits_read_key_lng(fits,(char*)name,&val,NULL, status)) {
        fits_report_error(stderr,*status);
    }
    return val;
}


int ft_get_colnum(fitsfile *fits, const char *colname, *status)
{
    int colnum=0;
    if (fits_get_colnum(fits, 0, (char*)colname, &colnum, status)) {
        fits_report_error(stderr,*status);
        colnum=FT_DEFVAL;
    }
    return colnum;
}

// size of a fixed-length array column (per row)
long ft_get_array_col_size(fitsfile *fits, const char *colname, int *status)
{
    int maxdim=2;
    int naxis = 0;
    long naxes[2]={0};

    int colnum=ft_get_colnum(fits, colname, status);
    if (*status != 0) {
        return FT_DEFVAL;
    }
    if (fits_read_tdim(fits, colnum, maxdim, &naxis,
                       naxes, status)) {
        fits_report_error(stderr,*status);
        return FT_DEFVAL;
    }
    return naxes[0];
}


char *ft_read_str_row_byname(fitsfile *fits,
                             const char* colname,
                             long row,
                             int *status)
{

    int colnum=ft_get_colnum(fits,colname,status);
    if (*status != 0) {
        (*status)=1;
        return NULL;
    }

    long len = ft_get_array_col_size(fits,colname);
    char* nulstr=" ";
    LONGLONG firstelem=1;

    char *data=calloc(len, sizeof(char));
    if (!data) {
        fprintf(stderr,"Could not allocate string\n");
        exit(1);
    }

    if (fits_read_col_str(fits, colnum, row, firstelem, 1,
                          nulstr, &data, NULL, status)) {
        fits_report_error(stderr,(*status));
        free(data);
    }
    return data;
}
static double fits_read_dbl_row_byname(fitsfile *fits,
                                       const char* colname, 
                                       long row,
                                       int *status)
{
    double nulval=0;
    LONGLONG firstelem=1;

    int colnum=get_colnum(fits,colname);

    if (*status != 0) {
        return FT_DEFVAL;
    }

    double data=0;
    if (fits_read_col_dbl(fits, colnum, row, firstelem, 1,
                          nulval, &data, NULL, status)) {
        fits_report_error(stderr,(*status));
        data=FT_DEFVAL;
    }
    return data;
}


