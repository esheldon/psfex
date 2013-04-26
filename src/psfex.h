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

#define PSFEX_NEIGEN(im) ((im)->neigen)

#define PSFEX_SIZE_PER(im) ((im)->eigen_size)
#define PSFEX_NROW_PER(im) ((im)->eigen_nrow)
#define PSFEX_NCOL_PER(im) ((im)->eigen_ncol)

#define PSFEX_SIZE_TOT(im) ((im)->mosaic_size)
#define PSFEX_NROW_TOT(im) ((im)->mosaic_nrow)
#define PSFEX_NCOL_TOT(im) ((im)->mosaic_ncol)

#define PSFEX_GET(im, eigen, row, col)                       \
    ( *((im)->rows[(eigen)*(im)->eigen_nrow + (row)] + (col)) )

#define PSFEX_GET_ROW(im, eigen, row)                       \
    ((im)->rows[(eigen)*(im)->eigen_nrow + (row)] )

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


struct psfex_eigens *psfex_eigens_new(long neigen,
                                      long nrow,
                                      long ncol);

struct psfex_eigens *psfex_eigens_free(struct psfex_eigens *self);
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
void psfex_write(const struct psfex *self, FILE* stream);

#endif
