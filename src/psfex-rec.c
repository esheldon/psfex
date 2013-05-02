#include <stdlib.h>
#include <stdio.h>
#include "psfex.h"
#include "psfex_fits.h"

int main(int argc, char **argv)
{
    if (argc < 5) {
        printf("usage: psfex-rec psfex_file outfile row col\n");
        printf("   reconstruct the psf image at the indicated location\n");
        printf("   and write to the output fits file\n");
        return 1;
    }

    const char *fname=argv[1];
    const char *outname=argv[2];
    double row = atof(argv[3]);
    double col = atof(argv[4]);

    int clobber=1;
    int status=0;

    struct psfex *psfex=psfex_fits_read(fname);

    struct psfex_image *im=psfex_rec_image(psfex, row, col);

    psfex_image_write_fits(im, outname, clobber, &status);

    im=psfex_image_free(im);
    psfex=psfex_free(psfex);
}
