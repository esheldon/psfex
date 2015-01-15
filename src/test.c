#include <stdlib.h>
#include <stdio.h>
#include "psfex.h"
#include "psfex_fits.h"

int main(int argc, char **argv)
{

    if (argc < 2) {
        printf("usage: test psfex_file [image_name]\n");
        printf("   optionally write out a reconstructed psf image\n");
        return 1;
    }

    const char *fname=argv[1];
    const char *outname=NULL;
    if (argc > 2) {
        outname=argv[2];
    }

    struct psfex *psfex=psfex_fits_read(fname);
    psfex_write(psfex, stdout);

    long row=11, col=8;
    long ncomp = PSFEX_NCOMP(psfex);
    for (long comp=0; comp < ncomp; comp++) {
	printf("Comp: %ld\n",comp);
	printf("  pix[%ld,%ld]: %lf\n", row, col, PSFEX_GET(psfex,comp,row,col));
    }

    double trow = 500.75;
    double tcol = 600.3;
    
    struct psfex_image *im=psfex_rec_image(psfex,trow, tcol);
    printf("rec[%ld,%ld]: %lf\n", row, col, PSFIM_GET(im, row, col));

    double rowcen,colcen;
    get_center(RECON_NROW(psfex),RECON_NCOL(psfex),trow,tcol,psfex->pixstep,&rowcen,&colcen);


    printf("center: row: %lf, col: %lf\n", rowcen, colcen);
 
    if (outname) {
        int clobber=1;
        int status=0;
        printf("writing fits: %s\n", outname);
        psfex_image_write_fits(im, outname, clobber, &status);
    }
    im=psfex_image_free(im);
    psfex=psfex_free(psfex);
}
