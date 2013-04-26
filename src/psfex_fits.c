#include <fitsio.h>
#include "psfex_fits.h"

static void read_keys(fitsfile *fits,
                      long *neigen,
                      long *nrow,   // per eigen
                      long *ncol,   // per eigen
                      long *poldeg,
                      double *polzero_row,
                      double *polzero_col,
                      double *polscale_row,
                      double *polscale_col,
                      double *psf_samp,
                      int *status)
{
    // first some checks
    long polnaxis;
    if (fits_read_key_lng(fits,"POLNAXIS",&polnaxis,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (polnaxis != 2) {
        fprintf(stderr,"expected POLNAXIS==2, got %ld", polnaxis);
        *status=1;
        return;
    }
    long psfnaxis;
    if (fits_read_key_lng(fits,"PSFNAXIS",&psfnaxis,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (psfnaxis != 3) {
        fprintf(stderr,"expected PSFNAXIS==3, got %ld", psfnaxis);
        *status=1;
        return;
    }
    long polngrp;
    if (fits_read_key_lng(fits,"POLNGRP",&polngrp,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (polngrp != 1) {
        fprintf(stderr,"expected POLNGRP==1, got %ld", polngrp);
        *status=1;
        return;
    }



    // these we keep
    if (fits_read_key_lng(fits,"PSFAXIS3",neigen,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (fits_read_key_lng(fits,"PSFAXIS2",nrow,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (fits_read_key_lng(fits,"PSFAXIS1",ncol,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (fits_read_key_lng(fits,"POLDEG1",poldeg,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }


    if (fits_read_key_dbl(fits,"POLZERO2",polzero_row,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (fits_read_key_dbl(fits,"POLZERO1",polzero_col,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (fits_read_key_dbl(fits,"POLSCAL2",polscale_row,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (fits_read_key_dbl(fits,"POLSCAL1",polscale_col,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }

    if (fits_read_key_dbl(fits,"PSF_SAMP",psf_samp,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }



}

static void read_eigens(struct psfex *self, fitsfile *fits, int *status)
{

    double nulval=0;
    LONGLONG firstrow=1;
    LONGLONG firstelem=1;
    LONGLONG nread=self->eigens->mosaic_size;

    double *data=self->eigens->rows[0];

    int colnum=0;
    if (fits_get_colnum(fits, 0, "PSF_MASK", &colnum, status)) {
        fits_report_error(stderr,*status);
        return;
    }
 
    if (fits_read_col_dbl(fits, colnum, firstrow, firstelem, nread,
                          nulval, data, NULL, status)) {
        fits_report_error(stderr,(*status));
    }
}
static struct psfex *psfex_from_fits(fitsfile *fits)
{
    struct psfex *self=NULL;
    long neigen;
    long nrow=0;   // per eigen
    long ncol=0;   // per eigen
    long poldeg=0;
    double polzero_row=0;
    double polzero_col=0;
    double polscale_row=0;
    double polscale_col=0;
    double psf_samp=0;
    int status=0;

    read_keys(fits,
              &neigen,
              &nrow,   // per eigen
              &ncol,   // per eigen
              &poldeg,
              &polzero_row,
              &polzero_col,
              &polscale_row,
              &polscale_col,
              &psf_samp,
              &status);
    if (status != 0 ) {
        goto _psfex_from_fits_bail;
    }

    self=psfex_new(neigen,
                   nrow,   // per eigen
                   ncol,   // per eigen
                   poldeg,
                   polzero_row,
                   polzero_col,
                   polscale_row,
                   polscale_col,
                   psf_samp);

    read_eigens(self, fits, &status);
    if (status != 0) {
        self=psfex_free(self);
    }

_psfex_from_fits_bail:
    return self;

}
struct psfex *psfex_fits_read(const char *filename)
{
    int status=0;
    fitsfile *fits=NULL;
    struct psfex *self=NULL;

    if (fits_open_file(&fits, filename, READONLY, &status)) {
        fits_report_error(stderr,status);
        goto _psex_fits_read_bail;
    }

    if (fits_movnam_hdu(fits, BINARY_TBL, "PSF_DATA", 0, &status)) {
        fits_report_error(stderr,status);
        goto _psex_fits_read_bail;
    }

    self = psfex_from_fits(fits);

_psex_fits_read_bail:
    if (fits) {
        if (fits_close_file(fits, &status)) {
            fits_report_error(stderr,status);
        }
    }
    return self;
}
