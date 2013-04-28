/*
   Based on PSFEx.h by peter melchior
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psfex.h"

static const double INTERPFAC = 3.0;
static const double IINTERPFAC = .3333333333333333333333333333;

static double sinc(double x) {
    if (x<1e-5 && x>-1e-5)
        return 1.;
    return sin(x*M_PI)/(x*M_PI);
}

static
struct psfex_eigens *psfex_eigens_new(long neigen,
                                      long nrow,   // per eigen
                                      long ncol)   // per eigen
{
    struct psfex_eigens *self=calloc(1, sizeof(struct psfex_eigens));
    if (!self) {
        fprintf(stderr,"failed to allocate struct psfex_eigens\n");
        exit(1);
    }

    self->neigen=neigen;

    self->mosaic_size = neigen*nrow*ncol;
    self->mosaic_nrow = neigen*nrow;
    self->mosaic_ncol = ncol;

    self->eigen_size = nrow*ncol;
    self->eigen_nrow=nrow;
    self->eigen_ncol=ncol;

    self->rows = calloc(self->mosaic_nrow,sizeof(double *));
    if (!self->rows) {
        fprintf(stderr,"could not allocate %ld image rows\n",
                self->mosaic_nrow);
        exit(1);
    }

    self->rows[0]=calloc(self->mosaic_nrow*self->mosaic_ncol,sizeof(double));

    for(long i = 1; i < self->mosaic_nrow; i++) {
        self->rows[i] = self->rows[i-1] + self->eigen_ncol;
    }

    return self;
}

static
struct psfex_eigens *psfex_eigens_free(struct psfex_eigens *self)
{
    if (self) {
        if (self->rows) {
            if (self->rows[0]) {
                free(self->rows[0]);
            }
            self->rows[0]=NULL;

            free(self->rows);
            self->rows=NULL;
        }
        free(self);
        self=NULL;
    }
    return self;
}

struct psfex *psfex_new(long neigen,
                        long nrow,   // per eigen
                        long ncol,   // per eigen
                        long poldeg,
                        double polzero_row,
                        double polzero_col,
                        double polscale_row,
                        double polscale_col,
                        double psf_samp)
{
    struct psfex *self=calloc(1, sizeof(struct psfex));
    if (!self) {
        fprintf(stderr,"failed to allocate struct psfex\n");
        exit(1);
    }

    self->poldeg=poldeg;
    self->polzero_row=polzero_row;
    self->polzero_col=polzero_col;
    self->polscale_row=polscale_row;
    self->polscale_col=polscale_col;
    self->psf_samp=psf_samp;

    self->eigens=psfex_eigens_new(neigen, nrow, ncol);
    if (!psfex_check(self)) {
        self=psfex_free(self);
        return self;
    }

    // maximum radius in the sample space (x/psf_samp)
    self->maxrad = (ncol-1)/2. - INTERPFAC;

    return self;
}

struct psfex *psfex_free(struct psfex *self)
{
    if (self) {
        if (self->eigens) {
            self->eigens=psfex_eigens_free(self->eigens);
        }
        free(self);
        self=NULL;
    }
    return self;
}

int psfex_check(const struct psfex *self)
{
    long neigen_exp = ((self->poldeg+1)*(self->poldeg+2))/2;
    if (self->eigens->neigen != neigen_exp) {
        fprintf(stderr,"poldeg and neigen disagree\n");
        return 0;
    }

    return 1;
}

void psfex_write(const struct psfex *self, FILE* stream)
{
    fprintf(stream,"poldeg:       %ld\n", self->poldeg);
    fprintf(stream,"polzero_row:  %lf\n", self->polzero_row);
    fprintf(stream,"polzero_col:  %lf\n", self->polzero_col);
    fprintf(stream,"polscale_row: %lf\n", self->polscale_row);
    fprintf(stream,"polscale_col: %lf\n", self->polscale_col);
    fprintf(stream,"psf_samp:     %lf\n", self->psf_samp);
    fprintf(stream,"maxrad:       %lf\n", self->maxrad);
    fprintf(stream,"neigen:       %ld\n", self->eigens->neigen);
    fprintf(stream,"nrow:         %ld\n", self->eigens->eigen_nrow);
    fprintf(stream,"ncol:         %ld\n", self->eigens->eigen_ncol);
}

struct psfex_image *psfex_image_new(long nrow, long ncol)
{
    return _psfex_image_new(nrow, ncol, 1);
}
struct psfex_image *_psfex_image_new(long nrow, long ncol, int alloc_data)
{
    struct psfex_image *self=calloc(1, sizeof(struct psfex_image));
    if (!self) {
        fprintf(stderr,"Could not allocate struct psfex_image\n");
        exit(1);
    }
    
    self->size=nrow*ncol;
    self->nrow=nrow;
    self->ncol=ncol;

    self->rows = calloc(self->nrow,sizeof(double *));
    if (!self->rows) {
        fprintf(stderr,"could not allocate %ld image rows\n",
                self->nrow);
        exit(1);
    }

    if (alloc_data) {
        self->rows[0]=calloc(self->size,sizeof(double));
        if (self->rows[0]==NULL) {
            fprintf(stderr,"could not allocate image of dimensions [%lu,%lu]\n",
                    nrow,ncol);
            exit(1);
        }

        for(long i = 1; i < self->nrow; i++) {
            self->rows[i] = self->rows[i-1] + self->ncol;
        }
        self->is_owner=1;
    } else {
        self->rows[0] = NULL;
        self->is_owner=0;
    }

    return self;
}
struct psfex_image *psfex_image_free(struct psfex_image *self)
{
    if (self) {
        if (self->rows) {
            if (self->rows[0]) {
                free(self->rows[0]);
            }
            self->rows[0]=NULL;

            free(self->rows);
            self->rows=NULL;
        }
        free(self);
        self=NULL;
    }
    return self;
}



/*

   Add up the contributions for each eigen image

   The row_scaled and col_scaled are the central row and col in the
   translated and scaled coordinates for the polynomial, (row-zero_row)/scale

   The erow, ecol are the pixel coords for the eigen images
*/
static
double get_summed_eigen_pixel(const struct psfex *self,
                              double row_scaled, double col_scaled, 
                              long erow, long ecol)
{
    // always start with value in the zeroth eigenimage
    double res=PSFEX_GET(self, 0, erow, ecol);

    for (long p=1; p<self->poldeg; p++) {
        for (long prow=0; prow<=p; prow++) {
            long pcol = p-prow;
            long k = pcol+prow*(self->poldeg+1)-(prow*(prow-1))/2;
            double eigval = PSFEX_GET(self, k, erow, ecol);
            res += pow(row_scaled,pcol) * pow(col_scaled,prow)* eigval;
        }
    }
    return res;
}
/*

   Get the pixel value.  The central row and col are in the translated and
   scaled coordinates for the polynomial, (row-zero_row)/scale

   The drow_samp,dcol_samp are relative to the *unscaled* centroid but are
   corrected for the sampling, e.g. (rowpsf-row)/psf_samp

   We then interpolate the pixels from a neighborhood radius defined (in the
   sample scale corrected coords) INTERPFAC

   you should check against maxrad before calling this function
*/
static
double get_pixel_value_samp(const struct psfex *self,
                            double row_scaled, double col_scaled,
                            double drow_samp, double dcol_samp)
{
    double pixval=0;

    long nrow=PSFEX_NROW(self);
    long ncol=PSFEX_NCOL(self);

    // interpolate values from the eigen images
    // we limit to the region defined by INTERPFAC
    // erow,ecol is for row in the eigen image set
    for(long erow=0; erow<nrow; erow++) {
        double derow = fabs(erow - 0.5*nrow - drow_samp);
        if (derow > INTERPFAC)
            continue;

        double derowdiv = derow*IINTERPFAC;

        for(long ecol=0; ecol<ncol; ecol++) {
            double decol = fabs(ecol - 0.5*ncol - dcol_samp);
            if (decol > INTERPFAC)
                continue;

            double decoldiv = decol*IINTERPFAC;

            double interpolant = 
                sinc(derow)*sinc(derowdiv)*sinc(decol)*sinc(decoldiv);

            double value = get_summed_eigen_pixel(self, 
                                                  row_scaled, col_scaled,
                                                  erow, ecol);
            pixval+= value*interpolant;
        }
    }

    return pixval;
}

double *psfex_recp(const struct psfex *self,
                   double row,
                   double col,
                   long *nrow,
                   long *ncol)
{

    (*nrow) = PSFEX_NROW(self);
    (*ncol) = PSFEX_NCOL(self);
    long npix=(*nrow)*(*ncol);
    double *data=calloc(npix, sizeof(double));
    if (!data) {
        fprintf(stderr,"could not allocate %ld doubles\n", npix);
        exit(1);
    }

    _psfex_rec_fill(self, row, col, data);

    return data;
}


struct psfex_image *psfex_rec_image(const struct psfex *self,
                                    double row,
                                    double col)
{
    long nrow=0, ncol=0;
    double *data=psfex_recp(self, row, col, &nrow, &ncol);

    // 0 means don't allocate the data
    struct psfex_image *im=_psfex_image_new(nrow, ncol, 0);

    im->rows[0] = data;

    for(long i = 1; i < im->nrow; i++) {
        im->rows[i] = im->rows[i-1] + im->ncol;
    }
    im->is_owner=1;
    return im;
}

void _psfex_rec_fill(const struct psfex *self,
                     double row,
                     double col,
                     double *data)
{

    long nrow = PSFEX_NROW(self);
    long ncol = PSFEX_NCOL(self);

    double row_scaled = (row-self->polzero_row)/self->polscale_row;
    double col_scaled = (col-self->polzero_col)/self->polscale_col;

    double sampfac = 1./(self->psf_samp*self->psf_samp);

    double rowpsf_cen=(nrow-1.)/2.;
    double colpsf_cen=(ncol-1.)/2.;

    for (long rowpsf=0; rowpsf<nrow; rowpsf++) {
        double drow_samp = (rowpsf-rowpsf_cen)/self->psf_samp;
        if (fabs(drow_samp) > self->maxrad)
            continue;

        for (long colpsf=0; colpsf<ncol; colpsf++) {

            double dcol_samp = (colpsf-colpsf_cen)/self->psf_samp;
            if (fabs(dcol_samp) > self->maxrad)
                continue;

            // pixle value in sample coords
            double pixval = get_pixel_value_samp(self, row_scaled, col_scaled, 
                                                 drow_samp, dcol_samp);
            // in pixel coords
            pixval *= sampfac;

            data[rowpsf*ncol + colpsf] = pixval;
        }
    }
}



