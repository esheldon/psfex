#include <string.h>
#include <fitsio.h>
#include "psfex_fits.h"

#define CFITSIO_MAX_ARRAY_DIMS 99

static void read_keys(fitsfile *fits,
		      long *masksize,  // [MASK_DIM]
		      long *poldeg,
		      double *contextoffset, // [POLY_DIM]
		      double *contextscale, // [POLY_DIM]
		      double *psf_samp,
		      int *status)
{
    long polnaxis, psfnaxis, polngrp;
    int i;
    char gstr[MAX_STR];

    // First part: check for dimensionality...
    //  (In future we can port over the rest of the sextractor code
    //   to read in different types of psf files if this becomes necessary.)
    
    
    // Read in the polynomial number of axes (require 2)
    if (fits_read_key_lng(fits,"POLNAXIS",&polnaxis,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (polnaxis != POLY_DIM) {
	fprintf(stderr,"Expected POLNAXIS==%d, got %ld\n",POLY_DIM, polnaxis);
	*status=1;
	return;
    }

    // Read in the psf number of axes (require 3)
    if (fits_read_key_lng(fits,"PSFNAXIS",&psfnaxis,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (psfnaxis != MASK_DIM) {
        fprintf(stderr,"expected PSFNAXIS==%d, got %ld", MASK_DIM, psfnaxis);
        *status=1;
        return;
    }

    // Read in number of groups (require 1)
    if (fits_read_key_lng(fits,"POLNGRP",&polngrp,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }
    if (polngrp != POLY_NGROUP) {
        fprintf(stderr,"expected POLNGRP==%d, got %ld", POLY_NGROUP, polngrp);
        *status=1;
        return;
    }

    // These next values we have to record...
    for (i=0;i<MASK_DIM;i++) {
	snprintf(gstr,MAX_STR,"PSFAXIS%1d", i+1);
	if (fits_read_key_lng(fits, gstr, &masksize[i], NULL, status)) {
	    fits_report_error(stderr,*status);
	    return;
	}
    }
    
    if (fits_read_key_lng(fits,"POLDEG1",poldeg,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }

    for (i=0;i<POLY_DIM;i++) {
	snprintf(gstr,MAX_STR,"POLZERO%1d", i+1);
	if (fits_read_key_dbl(fits,gstr,&contextoffset[i],NULL,status)) {
	    fits_report_error(stderr,*status);
	    return;
	}
	snprintf(gstr,MAX_STR,"POLSCAL%1d", i+1);
	if (fits_read_key_dbl(fits,gstr,&contextscale[i],NULL,status)) {
	    fits_report_error(stderr,*status);
	    return;
	}
    }
    
    if (fits_read_key_dbl(fits,"PSF_SAMP",psf_samp,NULL, status)) {
        fits_report_error(stderr,*status);
        return;
    }

}

static void read_psf_mask(struct psfex *self, fitsfile *fits, int *status)
{
    double nulval=0;
    LONGLONG firstrow=1;
    LONGLONG firstelem=1;
    LONGLONG nread;

    int colnum=0;

    nread = PSFEX_SIZE(self) * PSFEX_NCOMP(self);
    
    if (fits_get_colnum(fits, 0, "PSF_MASK", &colnum, status)) {
	fits_report_error(stderr,*status);
	return;
    }

    if (fits_read_col_dbl(fits, colnum, firstrow, firstelem, nread,
			  nulval, self->maskcomp, NULL, status)) {
	fits_report_error(stderr,*status);
    }
}

static struct psfex *psfex_from_fits(fitsfile *fits)
{
    struct psfex *self=NULL;
    long masksize[MASK_DIM],poldeg;
    double contextoffset[POLY_DIM], contextscale[POLY_DIM], psf_samp;
    int status=0;

    read_keys(fits,
	      masksize,
	      &poldeg,
	      contextoffset,
	      contextscale,
	      &psf_samp,
	      &status);

    if (status != 0 ) {
        goto _psfex_from_fits_bail;
    }

    self=psfex_new(masksize,
		   poldeg,
		   contextoffset,
		   contextscale,
		   psf_samp);

    read_psf_mask(self, fits, &status);

    if (status != 0) {
        self=psfex_free(self);
    }

	
_psfex_from_fits_bail:
    return self;


}

/*
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
*/

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
