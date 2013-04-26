import fitsio
import numpy

INTERPFAC = 3.0
EIGEN_DTYPE = 'f8'

class PSFEx(object):
    def __init__(self, filename):
        self._filename=filename
        self._load()

        self._check()

        self._set_maxrad()
        self._set_eigen_list()
        #self._set_grid()

    def rec(self, row, col):
        """
        Reconstruct the PSF image at the specified location

        psf is a double array
        """

        # always start with exactly first eigen image
        psf=self._eigen[0].copy()

        for i in xrange(1,self._neigen):
            pass

    def _set_grid(self):
        rgrid, cgrid=numpy.mgrid[0:self._nrows, self._ncols]
        self._rowcen=(self._nrows-1)/2
        self._colcen=(self._ncols-1)/2
        maxrad=self._maxrad
        self._wgood, = \
            numpy.where(  (numpy.abs(self._rgrid-self._rowcen) < maxrad)
                        & (numpy.abs(self._cgrid-self._colcen) < maxrad) )

    def _load(self):
        """
        Load the PSF information from the fits file
        """

        with fitsio.FITS(self._filename) as fits:
            # eigen means "eigen images"
            eigen_tot=fits['psf_data']['psf_mask'][:]
            h=fits['psf_data'].read_header()

        self._eigen_tot=eigen_tot.astype(EIGEN_DTYPE)

        self._neigen = self._eigen_tot.shape[0]
        self._nrows=self._eigen_tot.shape[1]
        self._ncols=self._eigen_tot.shape[2]


        self._hdata=h

        self._poldeg = h['poldeg1']
        self._polzero_row = h['polzero2']
        self._polzero_col = h['polzero1']

        self._polscale_row = h['polscal2']
        self._polscale_col = h['polscal1']

        self._psf_samp = h['psf_samp']

    def _set_eigen_list(self):
        elist=[]
        for i in xrange(self._neigen):
            elist.append(self._eigen_tot[i,:,:])

        self._eigen=elist

    def _set_maxrad(self):
        self._maxrad = ((self._ncols-1)/2.-INTERPFAC)*self._psf_samp

    def _check(self):
        h=self._hdata

        polngrp=h['polngrp']
        if polngrp != 1:
            raise ValueError("expected POLNGRP==1 but got %d" % polngrp)
        psfnaxis=h['psfnaxis']
        if psfnaxis != 3:
            raise ValueError("expected PSFNAXIS==3 but got %d" % psfnaxis)

        polnaxis=h['polnaxis']
        if polnaxis != 2:
            raise ValueError("expected POLNAXIS==2 but got %d" % polnaxis)

        neigen=self._neigen
        poldeg=self._poldeg
        if neigen != ((poldeg+1)*(poldeg+2))/2:
            raise ValueError("poldeg and neigen disagree")
