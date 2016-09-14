import fitsio
import numpy
from . import _psfex_pywrap

DEFAULT_EXT='psf_data'

DOUBLE_DTYPE = 'f8'
LONG_DTYPE = 'i8'
POLY_DIM = 2
MASK_DIM = 3
POLY_NGROUP = 1

class PSFExError(Exception):
    """
    Something fundamentally wrong with the PSFEx file
    """
    def __init__(self, value):
         self.value = value
    def __str__(self):
        return repr(self.value)


class PSFEx(dict):
    def __init__(self, *args, **kw):
        """
        Load data into a PSFEx object

        parameters
        ----------
        filename or FITS object: string, optional
            Send filename or open fitsio FITS object
        ext: string or int, optional
            The extension from which to read
        """

        self.load(*args, **kw)

    def get_rec(self, row, col):
        """
        Reconstruct the PSF image at the specified location

        psf is a double array
        """

        return self._psfex.rec(row, col)

    def get_center(self, row, col):
        """
        Get the "center" in the reconstruction for the
        input location
        """
        return self._psfex.center(row, col)

    def get_fwhm(self):
        """
        Get the fwhm in pixels from the header
        """
        return self._hdata['psf_fwhm']

    def get_sigma(self):
        """
        Get the "sigma" in pixels from the header
        """
        fac = 2.3548200450309493
        fwhm=self.get_fwhm()
        sigma = fwhm/fac
        return sigma

    def load(self, arg, **kw):
        """
        Load data into a PSFEx object

        parameters
        ----------
        filename or FITS object: string, optional
            Send filename or open fitsio FITS object
        ext: string or int, optional
            The extension from which to read
        """

        if isinstance(arg, fitsio.FITS):
            self._load_from_fits(arg, **kw)
        else:
            self['filename'] = arg
            with fitsio.FITS(arg) as fits:
                self._load_from_fits(fits, **kw)
 
    def _load_from_fits(self, fits, ext=DEFAULT_EXT):

        hdu = fits[ext]

        psf_mask=hdu['psf_mask'][:]

        psf_mask=numpy.array(psf_mask, dtype=DOUBLE_DTYPE)

        h=hdu.read_header()

        self._psf_mask=psf_mask

        if h['polnaxis'] != POLY_DIM:
            raise PSFExError("Expected POLNAXIS==%d, got %d" % (POLY_DIM, h['polnaxis']))

        if h['psfnaxis'] != MASK_DIM:
            raise PSFExError("Expected PSFNAXIS==%d, got %d" % (MASK_DIM, h['psfnaxis']))

        if h['polngrp'] != POLY_NGROUP:
            raise PSFExError("Expected POLNGRP==%d, got %d" % (POLY_NGROUP, h['polngrp']))

        self['poldeg'] = h['poldeg1']

        self['masksize'] = numpy.zeros(MASK_DIM, dtype=LONG_DTYPE)
        for i in range(MASK_DIM):
            self['masksize'][i] = h['psfaxis%1d' % (i+1)]
        
        self['contextoffset'] = numpy.zeros(POLY_DIM,dtype=DOUBLE_DTYPE)
        self['contextscale'] = numpy.zeros(POLY_DIM,dtype=DOUBLE_DTYPE)
        for i in range(POLY_DIM):
            self['contextoffset'][i] = h['polzero%1d' % (i+1)]
            self['contextscale'][i] = h['polscal%1d' % (i+1)]

        self['psf_samp'] = h['psf_samp']

        self._hdata = h

        self._psfex = _psfex_pywrap.PSFEx(self['masksize'],
                                          self['poldeg'],
                                          self['contextoffset'],
                                          self['contextscale'],
                                          self['psf_samp'],
                                          self._psf_mask)
 
       
