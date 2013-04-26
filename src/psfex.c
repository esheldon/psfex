#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psfex.h"

static const double INTERPFAC = 3.0;
static const double IINTERPFAC = .3333333333333333333333333333;

/*
static double sinc(double x) {
    if (x<1e-5 && x>-1e-5)
        return 1.;
    return sin(x*M_PI)/(x*M_PI);
}
*/

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

    self->maxrad = ((ncol-1)/2.-INTERPFAC)*self->psf_samp;

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

// image sampled psf pixel at position relx,rely relative to the psf center at centerx,centery
// brightest pixel is at relx=rely=0
/*
static double pixel_lanczos(const struct psfex *self,
                            double rowcen,
                            double colcen,
                            double row,
                            double col)
{
    double drowp = (row-rowcen)/self->psf_samp;
    double dcolp = (col-colcen)/self->psf_samp;

    double sum = 0.;

    for(int irow=0; irow<self->nrow; irow++) {
        double drow = fabs(irow - 0.5*self->nrow - drowp);
        if (drow > INTERPFAC)
            continue;

        double drowdiv = drow*IINTERPFAC;

        for(int icol=0; icol<self->ncol; icol++) {
            double dcol = fabs(icol - 0.5*self->ncol - dcolp);
            if (dcol > INTERPFAC)
                continue;

            double dcoldiv = dcol*IINTERPFAC;

            double interpolant = 
                sinc(drow)*sinc(drowdiv)*sinc(dcol)*sinc(dcoldiv);

            double value = 
                pixel_sampled_interp(self,irow,icol,rowcen,colcen);

            sum += value*interpolant;
        }
    }
    sum /= (psf_samp*psf_samp);
    return sum;
}
*/

