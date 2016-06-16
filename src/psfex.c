/*
   Based on PSFEx.h by peter melchior
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "psfex.h"
#include "poly.h"

static const double INTERPFAC = 4.0;
static const double IINTERPFAC = 0.25;

struct psfex *psfex_new(long *masksize, // [MASK_DIM]
			long poldeg,
			double *contextoffset, // [POLY_DIM]
			double *contextscale,  // [POLY_DIM]
			double psf_samp)
{
    struct psfex *self;
    int i;
    long deg[POLY_MAXDIM], group[POLY_MAXDIM];

    if ((self = calloc(1, sizeof(struct psfex))) == NULL) {
        fprintf(stderr,"Failed to allocate struct psfex\n");
        exit(1);
    }
    self->maskcomp = NULL;

    // copy masksize
    memcpy(self->masksize, masksize, MASK_DIM * sizeof(long));

    // copy contextoffset
    memcpy(self->contextoffset, contextoffset, POLY_DIM * sizeof(double));

    // copy contextscale
    memcpy(self->contextscale, contextscale, POLY_DIM * sizeof(double));

    self->masknpix = self->masksize[0] * self->masksize[1];

    // set up poly struct
    for (i=0;i<POLY_DIM;i++) {
        group[i] = 1;
    }

    deg[0] = poldeg;
    self->poly = poly_init(group, POLY_DIM, deg, POLY_NGROUP);

    self->pixstep = 1./psf_samp;

    // allocate memory for mask ...

    if ((self->maskcomp = (double *) calloc(self->masknpix * self->masksize[2], sizeof(double))) == NULL) {
        self=psfex_free(self);
        fprintf(stderr,"Failed to allocate maskcomp\n");
        exit(1);
    }

    // and set the reconstruction size, using psf sampling factor.  Make sure it's odd
    self->reconsize[0] = (long) ceil((float) self->masksize[0] * (float) psf_samp);
    if ((self->reconsize[0] % 2) == 0) { self->reconsize[0]++; }
    self->reconsize[1] = (long) ceil((float) self->masksize[1] * (float) psf_samp);
    if ((self->reconsize[1] % 2) == 0) { self->reconsize[1]++; } 
    

    return self;
}
			

struct psfex *psfex_free(struct psfex *self)
{
    if (self) {
        if (self->maskcomp) {
            free(self->maskcomp);
            self->maskcomp = NULL;
        }

        poly_end(self->poly);
        free(self);
    }
    return self;
}


void psfex_write(const struct psfex *self, FILE* stream)
{
    fprintf(stream,"masksize[0]:       %ld\n", self->masksize[0]);
    fprintf(stream,"masksize[1]:       %ld\n", self->masksize[1]);
    fprintf(stream,"masksize[2]:       %ld\n", self->masksize[2]);
    fprintf(stream,"reconsize[0]:      %ld\n", self->reconsize[0]);
    fprintf(stream,"reconsize[1]:      %ld\n", self->reconsize[1]);
    fprintf(stream,"contextoffset[0]:  %lf\n", self->contextoffset[0]);
    fprintf(stream,"contextoffset[1]:  %lf\n", self->contextoffset[1]);
    fprintf(stream,"contextscale[0]:   %lf\n", self->contextscale[0]);
    fprintf(stream,"contextscale[1]:   %lf\n", self->contextscale[1]);
    fprintf(stream,"pixstep:           %lf\n", self->pixstep);
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


void get_center(long nrow, long ncol,
		double row, double col,
		double pixstep,
		double *rowcen, double *colcen)
{
    double drow, dcol;

    // this follows the new calculation in _psfex_rec_fill, in turn
    //   following SExtractor
    
    dcol = col - (long) (col+0.5);
    drow = row - (long) (row+0.5);


    *colcen = (double) (ncol/2) + dcol;
    *rowcen = (double) (nrow/2) + drow;

}
		

/* repurposed from sextractor image.c */
/*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SExtractor is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SExtractor is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		19/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

static int _psfex_vignet_resample(double *pix1, int w1, int h1,
                                  double *pix2, int w2, int h2,
                                  double dx, double dy, double step2)
{
    double	*mask,*maskt, xc1,xc2,yc1,yc2, xs1,ys1, x1,y1, x,y, dxm,dym,
            val, norm,
            *pix12, *pixin,*pixin0, *pixout,*pixout0;
    int		i,j,k,n,t, *start,*startt, *nmask,*nmaskt,
            ixs2,iys2, ix2,iy2, dix2,diy2, nx2,ny2, iys1a, ny1, hmw,hmh,
            ix,iy, ix1,iy1;


    /* Initialize destination buffer to zero */
    memset(pix2, 0, w2*h2*sizeof(double));

    xc1 = (double)(w1/2);	/* Im1 center x-coord*/
    xc2 = (double)(w2/2);	/* Im2 center x-coord*/
    xs1 = xc1 + dx - xc2*step2;	/* Im1 start x-coord */

    if ((int)xs1 >= w1)
        return -1;
    ixs2 = 0;			/* Int part of Im2 start x-coord */
    if (xs1<0.0)
    {
        dix2 = (int)(1-xs1/step2);
        /*-- Simply leave here if the images do not overlap in x */
        if (dix2 >= w2)
            return -1;
        ixs2 += dix2;
        xs1 += dix2*step2;
    }
    nx2 = (int)((w1-1-xs1)/step2+1);/* nb of interpolated Im2 pixels along x */
    if (nx2>(ix2=w2-ixs2))
        nx2 = ix2;
    if (nx2<=0)
        return -1;
    yc1 = (double)(h1/2);	/* Im1 center y-coord */
    yc2 = (double)(h2/2);	/* Im2 center y-coord */
    ys1 = yc1 + dy - yc2*step2;	/* Im1 start y-coord */
    if ((int)ys1 >= h1)
        return -1;
    iys2 = 0;			/* Int part of Im2 start y-coord */
    if (ys1<0.0)
    {
        diy2 = (int)(1-ys1/step2);
        /*-- Simply leave here if the images do not overlap in y */
        if (diy2 >= h2)
            return -1;
        iys2 += diy2;
        ys1 += diy2*step2;
    }
    ny2 = (int)((h1-1-ys1)/step2+1);/* nb of interpolated Im2 pixels along y */
    if (ny2>(iy2=h2-iys2))
        ny2 = iy2;
    if (ny2<=0)
        return -1;

    /* Set the yrange for the x-resampling with some margin for interpolation */
    iys1a = (int)ys1;		/* Int part of Im1 start y-coord with margin */
    hmh = INTERPW/2 - 1;		/* Interpolant start */
    if (iys1a<0 || ((iys1a -= hmh)< 0))
        iys1a = 0;
    ny1 = (int)(ys1+ny2*step2)+INTERPW-hmh;	/* Interpolated Im1 y size */
    if (ny1>h1)					/* with margin */
        ny1 = h1;
    /* Express everything relative to the effective Im1 start (with margin) */
    ny1 -= iys1a;
    ys1 -= (double)iys1a;

    /* Allocate interpolant stuff for the x direction */
    if ((mask = (double *) malloc(sizeof(double) * nx2 * INTERPW)) == NULL) /* Interpolation masks */
        return -1;
    if ((nmask = (int *) malloc(sizeof(int) * nx2)) == NULL) /* Interpolation mask sizes */
        return -1;
    if ((start = (int *) malloc(sizeof(int) * nx2)) == NULL) /* Int part of Im1 conv starts */
        return -1;

    /* Compute the local interpolant and data starting points in x */
    hmw = INTERPW/2 - 1;
    x1 = xs1;
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    for (j=nx2; j--; x1+=step2)
    {
        ix = (ix1=(int)x1) - hmw;
        dxm = ix1 - x1 - hmw;	/* starting point in the interpolation func */
        if (ix < 0)
        {
            n = INTERPW+ix;
            dxm -= (double)ix;
            ix = 0;
        }
        else
            n = INTERPW;
        if (n>(t=w1-ix))
            n=t;
        *(startt++) = ix;
        *(nmaskt++) = n;
        norm = 0.0;
        for (x=dxm, i=n; i--; x+=1.0)
            norm += (*(maskt++) = INTERPF(x));
        norm = norm>0.0? 1.0/norm : 1.0;
        maskt -= n;
        for (i=n; i--;)
            *(maskt++) *= norm;
    }

    if ((pix12 = (double *) calloc(nx2*ny1, sizeof(double))) == NULL) { /* Intermediary frame-buffer */
        return -1;
    }

    /* Make the interpolation in x (this includes transposition) */
    pixin0 = pix1+iys1a*w1;
    pixout0 = pix12;
    for (k=ny1; k--; pixin0+=w1, pixout0++)
    {
        maskt = mask;
        nmaskt = nmask;
        startt = start;
        pixout = pixout0;
        for (j=nx2; j--; pixout+=ny1)
        {
            pixin = pixin0+*(startt++);
            val = 0.0; 
            for (i=*(nmaskt++); i--;)
                val += *(maskt++)**(pixin++);
            *pixout = val;
        }
    }

    /* Reallocate interpolant stuff for the y direction */
    if ((mask = (double *) realloc(mask, sizeof(double) * ny2 * INTERPW)) == NULL) { /* Interpolation masks */
        return -1;
    }
    if ((nmask = (int *) realloc(nmask, sizeof(int) * ny2)) == NULL) { /* Interpolation mask sizes */
        return -1;
    }
    if ((start = (int *) realloc(start, sizeof(int) * ny2)) == NULL) { /* Int part of Im1 conv starts */
        return -1;
    }

    /* Compute the local interpolant and data starting points in y */
    hmh = INTERPW/2 - 1;
    y1 = ys1;
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    for (j=ny2; j--; y1+=step2)
    {
        iy = (iy1=(int)y1) - hmh;
        dym = iy1 - y1 - hmh;	/* starting point in the interpolation func */
        if (iy < 0)
        {
            n = INTERPW+iy;
            dym -= (double)iy;
            iy = 0;
        }
        else
            n = INTERPW;
        if (n>(t=ny1-iy))
            n=t;
        *(startt++) = iy;
        *(nmaskt++) = n;
        norm = 0.0;
        for (y=dym, i=n; i--; y+=1.0)
            norm += (*(maskt++) = INTERPF(y));
        norm = norm>0.0? 1.0/norm : 1.0;
        maskt -= n;
        for (i=n; i--;)
            *(maskt++) *= norm;
    }

    /* Make the interpolation in y  and transpose once again */
    pixin0 = pix12;
    pixout0 = pix2+ixs2+iys2*w2;
    for (k=nx2; k--; pixin0+=ny1, pixout0++)
    {
        maskt = mask;
        nmaskt = nmask;
        startt = start;
        pixout = pixout0;
        for (j=ny2; j--; pixout+=w2)
        {
            pixin = pixin0+*(startt++);
            val = 0.0; 
            for (i=*(nmaskt++); i--;)
                val += *(maskt++)**(pixin++);
            *pixout = val;
        }
    }

    /* Free memory */
    free(pix12);
    free(mask);
    free(nmask);
    free(start);

    return 0;
}

void _psfex_rec_fill(const struct psfex *self,
                     double row,
                     double col,
                     double *data)
{
    double pos[POLY_DIM];
    double *basis=NULL, fac;
    double *ppc=NULL, *pl=NULL;
    int i,n,p;
    double *maskloc=NULL;
    double dcol,drow;
    //double sum;


    if ((maskloc = (double *) calloc(self->masknpix, sizeof(double))) == NULL) {
        fprintf(stderr,"Could not allocate maskloc\n");
        exit(1);
    }

    pos[0] = col;
    pos[1] = row;
    for (i=0;i<POLY_DIM;i++) {
        pos[i] = (pos[i] - self->contextoffset[i])/self->contextscale[i];
    }

    poly_func(self->poly, pos);

    basis = self->poly->basis;
    ppc = self->maskcomp;

    for (n=PSFEX_NCOMP(self); n--; ) {
        pl = maskloc;
        fac = *(basis++);
        for (p=PSFEX_SIZE(self); p--;)
            *(pl++) += fac**(ppc++);
    }

    // following sextractor where deltax = mx - ix
    //                            mx is the center
    //                            ix is the integer (floor) center of the stamp
    // this does not follow sextractor which has an extra -1 because YOLO
    //   (Note that this shifts the cutout in the postage stamp; get_center()
    //   tells you where the center actually is.  Also note that using the
    //   -1 gives boundary problems on the reconstructed postage stamp,
    //   which is why it has been removed for 0.3.1 -- ESR

    dcol = col - (long) (col+0.5);
    drow = row - (long) (row+0.5);
    
    _psfex_vignet_resample(maskloc,
                           self->masksize[0],
                           self->masksize[1],
                           data,
                           self->reconsize[0],
                           self->reconsize[1],
                           -dcol*self->pixstep,
                           -drow*self->pixstep,
                           self->pixstep);
    
    // NOTE: this is not normalized to match SExtractor fits at the moment...
    // This will be updated when/if SExtractor is updated...

    free(maskloc);
}


double *psfex_recp(const struct psfex *self,
                   double row,
                   double col,
                   long *nrow,
                   long *ncol)
{

    // this is the size of the reconstructed image
    (*nrow) = RECON_NROW(self);
    (*ncol) = RECON_NCOL(self);
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



