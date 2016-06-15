psfex
=====

python and C codes to interpolate and reconstruct psfex images

using the python code
---------------------

```python
import psfex

row=514.25
col=610.00
pex = psfex.PSFEx(filename)
image = pex.get_rec(row, column)
```

using the C library
-------------------

```C
// this version uses the fits reading
#include "psfex.h"
#include "psfex_fits.h"

double row=514.25, col=610.;

struct psfex *psfex=psfex_fits_read(fname);
struct psfex_image *im=psfex_rec_image(psfex,row,col);

printf("rec[%ld,%ld]: %lf\n", row, col, PSFIM_GET(im, row, col));

im=psfex_image_free(im);
psfex=psfex_free(psfex);
```

using the standalone psfex-rec code
----------------------------------

This code will write out a fits file with the PSF image in it,
to be read by your favorite code.

```bash
psfex-rec psfex_file output_file row col
```

The image will be in output_file

installation of python code
----------------------------

```bash
git clone https://github.com/esheldon/psfex.git

cd psfex

# to install globally
python setup.py install

# to install in a particular place
python setup.py install --prefix=/some/path

# if you install in a prefix, make sure you
# add the /some/path/lib/python2.7/site-packages
# directory to your PYTHONPATH (replace python2.7
# with your python version)
```

installation of C library and standalone psfex-rec code
------------------------------------------------------

```bash
git clone https://github.com/esheldon/psfex.git

cd psfex

# to install globally
make install

# to install in a particular place
make install prefix=/some/path

# if you install in a prefix, make sure you
# add the /some/path/bin to your PATH.
#
# also /some/path/lib # directory to your LD_LIBRARY_PATH and
# LIBRARY_PATH.  Similarly add /some/path/include
# to C_LIBRARY_PATH and CPATH
```

dependencies for python library
-------------------------------

- numpy
- fitsio - https://github.com/esheldon/fitsio

dependencies for C library
-------------------------------

- cfitsio
