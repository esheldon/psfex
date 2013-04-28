#ifndef _PSFEX_FITS_HEADER_GUARD
#define _PSFEX_FITS_HEADER_GUARD

#include "psfex.h"

struct psfex *psfex_fits_read(const char *filename);

void psfex_image_write_fits(const struct psfex_image *self,
                            const char *filename,
                            int clobber,
                            int *status);

struct psfex_image *psfex_image_read_fits(const char *fname, int ext,
                                          int *status);
#endif
