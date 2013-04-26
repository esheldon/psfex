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

    psfex=psfex_free(psfex);
}
