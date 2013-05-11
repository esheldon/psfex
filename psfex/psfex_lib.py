import fitsio
import numpy
from . import _psfex_pywrap

INTERPFAC = 3.0
EIGEN_DTYPE = 'f8'

class PSFEx(dict):
    def __init__(self, filename):
        self._load(filename)

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
        from math import floor
        rowcen_int=(self['nrow']-1)//2
        colcen_int=(self['ncol']-1)//2

        row_remain=row-floor(row)
        col_remain=col-floor(col)

        rowcen = float(rowcen_int) + row_remain
        colcen = float(colcen_int) + col_remain
        return rowcen, colcen

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

