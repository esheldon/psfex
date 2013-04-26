#include <stdlib.h>
#include <stdio.h>
#include "psfex.h"
#include "psfex_fits.h"

int main(int argc, char **argv)
{

    if (argc < 2) {
        printf("usage: test psfex_file\n");
        return 1;
    }

    const char *fname=argv[1];

    struct psfex *psfex=psfex_fits_read(fname);
    psfex_write(psfex, stdout);

    long row=11, col=8;
    long neigen=PSFEX_NEIGEN(psfex);
    for (long eigen=0; eigen<neigen; eigen++) {
        printf("eigen: %ld\n", eigen);
        printf("  pix[%ld,%ld]: %lf\n", row, col, PSFEX_GET(psfex,eigen,row,col));
    }

    psfex=psfex_free(psfex);
}
