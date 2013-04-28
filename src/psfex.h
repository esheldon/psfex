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

    struct psfex_eigens *eigens;
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

#define PSFIM_SIZE(im) ((im)->size)
#define PSFIM_NROW(im) ((im)->nrow)
#define PSFIM_NCOL(im) ((im)->ncol)
#define PSFIM_GET(im, row, col)      \
    ( *((im)->rows[(row)] + (col)) )

#define PSFIM_GETP(im, row, col)                 \
    (  ((im)->rows[(row)] + (col)) )

struct psfex_image *psfex_image_new(long nrow, long ncol);
struct psfex_image *_psfex_image_new(long nrow, long ncol, int alloc_data);

struct psfex_image *psfex_image_free(struct psfex_image *self);

// reconstruct a psf image at the indicated location
// and return a psfex_image object
struct psfex_image *psfex_rec_image(const struct psfex *self,
                                    double row,
                                    double col);



#endif
