#include <string.h>
#include <fitsio.h>
#include "psfex_fits.h"

#define CFITSIO_MAX_ARRAY_DIMS 99

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





/* 
   add a ! to front of name so cfitsio will clobber any existing file 
   you must free the returned string.
*/
static char *get_clobber_name(const char *filename)
{
    char *oname=NULL;
    int len=strlen(filename);

    oname = calloc(len+2, sizeof(char));
    oname[0]='!';

    strncpy(oname+1, filename, len);
    return oname;
}


void psfex_image_write_fits(const struct psfex_image *self,
                            const char *filename,
                            int clobber,
                            int *status)
{
    fitsfile* fits=NULL;
    LONGLONG firstpixel=1;
    LONGLONG nelements=0;

    int ndims=2;
    long dims[2]={0};

    char *name=NULL;

    if (clobber) {
        name=get_clobber_name(filename);
    } else {
        name=strdup(filename);
    }

    if (fits_create_file(&fits, name, status)) {
        fits_report_error(stderr,*status);
        goto _psfex_image_write_fits_bail;
    }

    dims[1] = PSFIM_NROW(self);
    dims[0] = PSFIM_NCOL(self);
    if (fits_create_img(fits, DOUBLE_IMG, ndims, dims, status)) {
        fits_report_error(stderr,*status);
        goto _psfex_image_write_fits_bail;
    }

    nelements=PSFIM_SIZE(self);
    if (fits_write_img(fits, TDOUBLE, firstpixel, nelements, 
                       PSFIM_GETP(self,0,0), status)) {
        fits_report_error(stderr,*status);
        goto _psfex_image_write_fits_bail;
    }

    if (fits_close_file(fits, status)) {
        fits_report_error(stderr,*status);
        goto _psfex_image_write_fits_bail;
    }

_psfex_image_write_fits_bail:
    free(name);
}

struct psfex_image *psfex_image_read_fits(const char *fname, int ext,
                                          int *status)
{
    fitsfile* fits=NULL;
    struct psfex_image *image=NULL;

    if (fits_open_file(&fits, fname, READONLY, status)) {
        fits_report_error(stderr, *status);
        goto _psfex_image_read_fits_bail;
    }

    int hdutype=0;
    if (fits_movabs_hdu(fits, ext+1, &hdutype, status)) {
        fits_report_error(stderr, *status);
        goto _psfex_image_read_fits_bail;
    }


    int maxdim=CFITSIO_MAX_ARRAY_DIMS;
    LONGLONG dims[CFITSIO_MAX_ARRAY_DIMS];
    int bitpix=0, ndims=0;
    if (fits_get_img_paramll(fits, maxdim, &bitpix, &ndims, dims, status)) {
        fits_report_error(stderr, *status);
        goto _psfex_image_read_fits_bail;
    }
    if (ndims != 2) {
        fprintf(stderr,"expected ndims=2, got %d\n", ndims);
        goto _psfex_image_read_fits_bail;
    }
    // dims reversed
    //fprintf(stderr,"dims: [%lld,%lld]\n", dims[0], dims[1]);

    // note dims are reversed
    image=psfex_image_new(dims[1], dims[0]);
    long npix=dims[1]*dims[0];
    long fpixel[2]={1,1};
    if (fits_read_pix(fits, TDOUBLE, fpixel, npix,
                      NULL,image->rows[0], NULL, status)) {
        fits_report_error(stderr, *status);
        goto _psfex_image_read_fits_bail;
    }

_psfex_image_read_fits_bail:
    if (*status) {
        image=psfex_image_free(image);
    }
    if (fits_close_file(fits, status)) {
        fits_report_error(stderr, *status);
    }

    return image;
}
