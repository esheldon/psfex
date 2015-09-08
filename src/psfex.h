#ifndef _PSFEX_HEADER_GUARD
#define _PSFEX_HEADER_GUARD


#define POLY_DIM    2
#define MASK_DIM    3
#define POLY_NGROUP 1
#define RECON_DIM   2

#define MAX_STR     100


struct psfex {
    // PSFEx/SExtractor based ...
    long         masksize[MASK_DIM];          /* Size of each dimension: fixed at 3 */
    long         masknpix;             /* Total number of PSF pixels */
    double      *maskcomp;             /* Complete pixel data (PSF components) */
    double       contextoffset[POLY_DIM];     /* Offset to apply to context data */
    double       contextscale[POLY_DIM];      /* Scale to apply to context data */
    struct poly *poly;                 /* Polynomial structure */
    double       pixstep;              /* PSF Sampling step */
    long         reconsize[RECON_DIM];         /* size of reconstructed image (after sampling) */
};

#define PSFEX_NCOMP(im) ((im)->masksize[2])

#define PSFEX_SIZE(im) ((im)->masknpix)
#define PSFEX_NROW(im) ((im)->masksize[1])
#define PSFEX_NCOL(im) ((im)->masksize[0])
#define RECON_NROW(im) ((im)->reconsize[1])
#define RECON_NCOL(im) ((im)->reconsize[0])

#define PSFEX_GET(im, comp, row, col)               \
    ( (im)->maskcomp[comp*PSFEX_NROW(im)*PSFEX_NCOL(im) + row*PSFEX_NCOL(im) + col])

struct psfex *psfex_new(long *masksize,
			long poldeg,
			double *contextoffset,
			double *contextscale,
			double psf_samp);
			

struct psfex *psfex_free(struct psfex *self);

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

//void get_center(long nrow, long ncol,
//		double row, double col,
//		double *rowcen, double *colcen);

void get_center(long nrow, long ncol,
		double row, double col,
		double pixstep,
		double *rowcen, double *colcen);


#endif
