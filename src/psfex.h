#ifndef _PSFEX_HEADER_GUARD
#define _PSFEX_HEADER_GUARD



struct psfex_eigens {
    long neigen;

    long mosaic_size;       // neigen*nrow*ncol
    long mosaic_nrow;       // neigen*nrow
    long mosaic_ncol;       // ncol

    long eigen_size;       // nrow*ncol
    long eigen_nrow;       // per eigen
    long eigen_ncol;       // per eigen

    double **rows;
};

#define PSFEX_NEIGEN(im) ((im)->eigens->neigen)

#define PSFEX_SIZE(im) ((im)->eigens->eigen_size)
#define PSFEX_NROW(im) ((im)->eigens->eigen_nrow)
#define PSFEX_NCOL(im) ((im)->eigens->eigen_ncol)

#define PSFEX_SIZE_TOT(im) ((im)->eigens->mosaic_size)
#define PSFEX_NROW_TOT(im) ((im)->eigens->mosaic_nrow)
#define PSFEX_NCOL_TOT(im) ((im)->eigens->mosaic_ncol)

#define PSFEX_GET(im, eigen, row, col)                       \
    ( *((im)->eigens->rows[(eigen)*(im)->eigens->eigen_nrow + (row)] + (col)) )

#define PSFEX_GET_ROW(im, eigen, row)                       \
    ((im)->eigens->rows[(eigen)*(im)->eigens->eigen_nrow + (row)] )

struct psfex {
    long poldeg;

    double polzero_row;
    double polzero_col;
    double polscale_row;
    double polscale_col;
    double psf_samp;

    double maxrad;

    double interpfac;
    double iinterpfac;

    struct psfex_eigens *eigens;

  // and the alternative stuff from sextractor...
  int		maskdim;	/* Dimensionality of the tabulated data */
  //int		*masksize;	/* PSF mask dimensions */
  int           masksize[2];
  int		masknpix;	/* Total number of involved PSF pixels */
  float		*maskcomp;      /* Complete pix. data (PSF components) */
  float		*maskloc;	/* Local PSF */
  //  double	**context;	/* Contexts */
  // t_type	*contexttyp;	/* Context types */
  // char		**contextname;	/* Array of context key-names */
  //double	*contextoffset;	/* Offset to apply to context data */
  //double	*contextscale;	/* Scaling to apply to context data */
  double        contextoffset[2];
  double        contextscale[2];
  struct poly	*poly;		/* Polynom describing the PSF variations */
  //pcstruct	*pc;		/* PC components */
  //double	fwhm;		/* Typical PSF FWHM */
  float		pixstep;	/* PSF sampling step */
  int		build_flag;	/* Set if the current PSF has been computed */
  
};


struct psfex *psfex_new(long neigen,
                        long nrow,   // per eigen
                        long ncol,   // per eigen
                        long poldeg,
                        double polzero_row,
                        double polzero_col,
                        double polscale_row,
                        double polscale_col,
                        double psf_samp);

struct psfex *psfex_free(struct psfex *self);
int psfex_check(const struct psfex *self);

// write the metadata to the input stream
void psfex_write(const struct psfex *self, FILE* stream);

// reconstruct a psf image at the indicated location
// return a pointer to allocated data
double *psfex_recp(const struct psfex *self,
                   double row,
                   double col,
                   long *nrow,
                   long *ncol);

// the user will probably want to use their own image class
// this is available if wanted and for testing purposes
struct psfex_image {
    size_t size;   // masked size
    size_t nrow;  // masked nrows
    size_t ncol;  // masked ncols

    int is_owner;
    double **rows;
};

// fill the user-supplied data.  The user is responsible
// for making sure the data is nrow*ncol in size!
void _psfex_rec_fill(const struct psfex *self,
                     double row,
                     double col,
                     double *data);

#define PSFIM_SIZE(im) ((im)->size)
#define PSFIM_NROW(im) ((im)->nrow)
#define PSFIM_NCOL(im) ((im)->ncol)
#define PSFIM_GET(im, row, col)      \
    ( *((im)->rows[(row)] + (col)) )

#define PSFIM_GETP(im, row, col)                 \
    (  ((im)->rows[(row)] + (col)) )

#define         INTERPW		8	/* Interpolation function range */
#define PI      	3.1415926535898

#define	INTERPF(x)	(x<1e-5 && x>-1e-5? 1.0 \
			:(x>INTERPFAC?0.0:(x<-INTERPFAC?0.0 \
			:sinf(PI*x)*sinf(PI/INTERPFAC*x)/(PI*PI/INTERPFAC*x*x))))
				/* Lanczos approximation */

struct psfex_image *psfex_image_new(long nrow, long ncol);
struct psfex_image *_psfex_image_new(long nrow, long ncol, int alloc_data);

struct psfex_image *psfex_image_free(struct psfex_image *self);

// reconstruct a psf image at the indicated location
// and return a psfex_image object
struct psfex_image *psfex_rec_image(const struct psfex *self,
                                    double row,
                                    double col);

int _psfex_vignet_resample(float *pix1, int w1, int h1,
			   float *pix2, int w2, int h2,
			   float dx, float dy, float step2);


#endif
