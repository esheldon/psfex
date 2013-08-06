/*
   Based on PSFEx.h by peter melchior
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "psfex.h"
#include "poly.h"

static const double INTERPFAC = 3.0;
static const double IINTERPFAC = .3333333333333333333333333333;

/*
static double sinc(double x) {
    if (x<1e-5 && x>-1e-5)
        return 1.;
    return sin(x*M_PI)/(x*M_PI);
}
*/

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

    self->rows[0]=calloc(self->mosaic_size,sizeof(double));

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
  int ndim, ngroup, deg[POLY_MAXDIM], group[POLY_MAXDIM];
  int psfnaxis;

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

    // and the sextractor stuff...
    ndim = 2;  // this is fixed
    ngroup = 1; // fixed
    group[0] = 1; // fixed
    group[1] = 1;
    psfnaxis = 3; // fixed

    self->contextoffset[0] = polzero_col;
    self->contextoffset[1] = polzero_row;
    self->contextscale[0] = polscale_col;
    self->contextscale[1] = polscale_row;
    self->masksize[0] = ncol;
    self->masksize[1] = nrow;
    self->masksize[2] = neigen;

    self->maskdim = psfnaxis;

    deg[0] = poldeg;
    //    fprintf(stdout,"About to do poly_init...\n");
    self->poly = poly_init(group,ndim,deg,ngroup);

    self->pixstep = (float) psf_samp;

    // read in eigens later...but allocate memory now.
    if ((self->maskcomp = (float *) calloc(neigen*nrow*ncol,sizeof(float))) == NULL) {
      self=psfex_free(self);
      return self;
    }

    // and allocated memory for maskloc
    if ((self->maskloc = (float *) calloc(nrow*ncol,sizeof(float))) == NULL) {
      self=psfex_free(self);
      return self;
    }

    //fprintf(stdout,"Done making new psf\n");

    return self;
}

struct psfex *psfex_free(struct psfex *self)
{
  //fprintf(stdout,"Freeing self\n");
    if (self) {
        if (self->eigens) {
            self->eigens=psfex_eigens_free(self->eigens);
        }
	poly_end(self->poly);
	free(self->maskcomp);
	free(self->maskloc);
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
/*
static
double get_summed_eigen_pixel(const struct psfex *self,
                              double row_scaled, double col_scaled, 
                              long erow, long ecol)
{
    // always start with value in the zeroth eigenimage
    double res=PSFEX_GET(self, 0, erow, ecol);

    for (long p=1; p<=self->poldeg; p++) {
        for (long prow=0; prow<=p; prow++) {
            long pcol = p-prow;
            long k = pcol+prow*(self->poldeg+1)-(prow*(prow-1))/2;
            double eigval = PSFEX_GET(self, k, erow, ecol);
            res += pow(col_scaled,pcol) * pow(row_scaled,prow)* eigval;
        }
    }
    return res;
}
*/
/*

   Get the pixel value.  The central row and col are in the translated and
   scaled coordinates for the polynomial, (row-zero_row)/scale

   The drow_samp,dcol_samp are relative to the *unscaled* centroid but are
   corrected for the sampling, e.g. (rowpsf-row)/psf_samp

   We then interpolate the pixels from a neighborhood radius defined (in the
   sample scale corrected coords) INTERPFAC

   you should check against maxrad before calling this function
*/
/*
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
*/

static void get_center(long nrow, long ncol,
                       double row, double col,
                       double *rowcen, double *colcen)
{

    long rowcen_int=(nrow-1)/2;
    long colcen_int=(ncol-1)/2;

    double row_remain=row-floor(row);
    double col_remain=col-floor(col);

    (*rowcen) = (double)rowcen_int + row_remain + 0.5;
    (*colcen) = (double)colcen_int + col_remain + 0.5;

}
void _psfex_rec_fill(const struct psfex *self,
                     double row,
                     double col,
                     double *data)
{
  // THIS IS THE MAIN CODE...

  //fprintf(stdout,"In _psfex_rec_fill\n");

  // replace with psf_build (!)
  // and fill the data array...this may be straightforward.

  static double pos[POLY_MAXDIM];
  double *basis, fac;
  float *ppc, *pl;
  int   i,n,p,ndim,npix;
  float *resampled;
  double rowpsf_cen=0, colpsf_cen=0;
  //int j,k;

  npix = self->masksize[0]*self->masksize[1];

  //fprintf(stdout,"memsetting...\n");
  memset(self->maskloc, 0, npix*sizeof(float));

  //fprintf(stdout,"getting pos...\n");
  ndim = self->poly->ndim;
  pos[0] = col;
  pos[1] = row;
  for (i=0;i<ndim;i++) {
    pos[i] = (pos[i] - self->contextoffset[i]) / self->contextscale[i];
  }

  //fprintf(stdout,"Running poly_func (%.2f/%.2f)...\n",pos[0],pos[1]);
  poly_func(self->poly, pos);

  /*
  fprintf(stdout,"POLY\nncoeff = %d\nndim = %d\nngroup=%d\n",self->poly->ncoeff,self->poly->ndim,self->poly->ngroup);
  fprintf(stdout,"BASIS:\n");
  for (i=0;i<self->poly->ndim;i++) {
    fprintf(stdout,"%f ",self->poly->basis[i]);
  }
  fprintf(stdout,"\nCOEFF:\n");
  for (i=0;i<self->poly->ndim;i++) {
    fprintf(stdout,"%f ",self->poly->coeff[i]);
  }
  fprintf(stdout,"\nDEGREE:\n");
  for (i=0;i<self->poly->ngroup;i++) {
    fprintf(stdout,"%d ",self->poly->degree[i]);
  }
  fprintf(stdout,"\n");
  */
  basis = self->poly->basis;

  ppc = self->maskcomp;

  //fprintf(stdout,"Craziness...\n");
  //fprintf(stdout,"maskdim, masksize[2] = %d, %d\n",self->maskdim, self->masksize[2]);
  for (n = (self->maskdim>2?self->masksize[2]:1); n--; ) {
    pl = self->maskloc;
    fac = *(basis++);
    //fprintf(stdout,"pl = %.5f, fac = %.5f, ppc = %.5f\n", pl[311], fac, ppc[311]);
    for (p=npix; p--;)
      *(pl++) += fac**(ppc++);
  }
  /*
  fprintf(stdout,"pos: %.2f, %.2f\n",pos[0],pos[1]);
  for (j=0;j<self->masksize[0];j++) {
    for (k=0;k<self->masksize[1];k++) {
      fprintf(stdout,"%.5f ",self->maskloc[k+j*self->masksize[0]]);
    }
    fprintf(stdout,"\n");
  }
  */

  //fprintf(stdout,"Copying data...\n");
  // and copy into the data...

  // Need to resample...
  if ((resampled = (float *) calloc(npix,sizeof(float))) == NULL) {
    return;
  }

  get_center(PSFEX_NROW(self), PSFEX_NCOL(self), row, col, &rowpsf_cen, &colpsf_cen);
  rowpsf_cen -= (float)PSFEX_NROW(self)/2;
  colpsf_cen -= (float)PSFEX_NCOL(self)/2;
  //fprintf(stdout,"center = %f, %f\n", rowpsf_cen, colpsf_cen);
  
  _psfex_vignet_resample(self->maskloc, self->masksize[0], self->masksize[1],
			 resampled, self->masksize[0], self->masksize[1],
			 -(float)colpsf_cen*self->pixstep, -(float)rowpsf_cen*self->pixstep,
			 self->pixstep);

  /*for (j=0;j<self->masksize[0];j++) {
    for (k=0;k<self->masksize[1];k++) {
      fprintf(stdout,"%.5f ",fabs(resampled[k+j*self->masksize[0]]));
    }
    fprintf(stdout,"\n");
    }*/

    
  for (i=0;i<npix;i++) {
    //  data[i] = (double) self->maskloc[i];
    data[i] = (double) resampled[i];
  }

  free(resampled);


  /*
    long nrow = PSFEX_NROW(self);
    long ncol = PSFEX_NCOL(self);

    double row_scaled = (row-self->polzero_row)/self->polscale_row;
    double col_scaled = (col-self->polzero_col)/self->polscale_col;

    double sampfac = 1./(self->psf_samp*self->psf_samp);

    double rowpsf_cen=0, colpsf_cen=0;
    get_center(nrow, ncol, row, col, &rowpsf_cen, &colpsf_cen);

    double sum=0;
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
            sum += pixval;
        }
    }

    double normfac=1./sum;
    long totpix=PSFEX_SIZE(self);
    for (long i=0; i<totpix; i++) {
        data[i] *= normfac;
    }
  */
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


/* from sextractor image.c */

int _psfex_vignet_resample(float *pix1, int w1, int h1,
			   float *pix2, int w2, int h2,
			   float dx, float dy, float step2)
{
  float	*mask,*maskt, xc1,xc2,yc1,yc2, xs1,ys1, x1,y1, x,y, dxm,dym,
    val, norm,
    *pix12, *pixin,*pixin0, *pixout,*pixout0;
  int		i,j,k,n,t, *start,*startt, *nmask,*nmaskt,
    ixs2,iys2, ix2,iy2, dix2,diy2, nx2,ny2, iys1a, ny1, hmw,hmh,
    ix,iy, ix1,iy1;
  

  /* Initialize destination buffer to zero */
  memset(pix2, 0, w2*h2*sizeof(float));
  
  xc1 = (float)(w1/2);	/* Im1 center x-coord*/
  xc2 = (float)(w2/2);	/* Im2 center x-coord*/
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
  yc1 = (float)(h1/2);	/* Im1 center y-coord */
  yc2 = (float)(h2/2);	/* Im2 center y-coord */
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
  ys1 -= (float)iys1a;

  /* Allocate interpolant stuff for the x direction */
  if ((mask = (float *) malloc(sizeof(float) * nx2 * INTERPW)) == NULL) /* Interpolation masks */
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
      dxm -= (float)ix;
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

  if ((pix12 = (float *) calloc(nx2*ny1, sizeof(float))) == NULL) { /* Intermediary frame-buffer */
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
  if ((mask = (float *) realloc(mask, sizeof(float) * ny2 * INTERPW)) == NULL) { /* Interpolation masks */
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
      dym -= (float)iy;
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
