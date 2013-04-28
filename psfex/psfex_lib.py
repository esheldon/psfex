import fitsio
import numpy
from . import _psfex_pywrap

INTERPFAC = 3.0
EIGEN_DTYPE = 'f8'

class PSFEx(dict):
    def __init__(self, filename):
        self._load(filename)

    def rec(self, row, col):
        """
        Reconstruct the PSF image at the specified location

        psf is a double array
        """

        return self._psfex.rec(row, col)

    def _load(self, filename):
        """
        Load the PSF information from the fits file

        Also load the C PSFEx object
        """

        self['filename'] = filename
        with fitsio.FITS(filename) as fits:
            eigens=fits['psf_data']['psf_mask'][:]
            eigens=numpy.array(eigens, dtype=EIGEN_DTYPE)

            h=fits['psf_data'].read_header()

        # make a copy, will be native byte order
        self._eigens=eigens

        self['neigen'] = self._eigens.shape[0]
        self['nrow']=self._eigens.shape[1]
        self['ncol']=self._eigens.shape[2]

        self._hdata=h

        self['poldeg'] = h['poldeg1']

        self['polzero_row'] = h['polzero2']
        self['polzero_col'] = h['polzero1']

        self['polscale_row'] = h['polscal2']
        self['polscale_col'] = h['polscal1']

        """
        self['polzero_row'] = h['polzero1']
        self['polzero_col'] = h['polzero2']

        self['polscale_row'] = h['polscal1']
        self['polscale_col'] = h['polscal2']
        """

        self['psf_samp'] = h['psf_samp']

        self._psfex = _psfex_pywrap.PSFEx(self['neigen'],
                                          self['nrow'],
                                          self['ncol'],
                                          self['poldeg'],
                                          self['polzero_row'],
                                          self['polzero_col'],
                                          self['polscale_row'],
                                          self['polscale_col'],
                                          self['psf_samp'],
                                          self._eigens)

