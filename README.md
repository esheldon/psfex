psfex
=====

python and C codes to interpolate and reconstruct psfex images

installation of python code
----------------------------

    git clone git@github.com:esheldon/psfex.git

    cd psfex

    # to install globally
    python setup.py install

    # to install in a particular place
    python setup.py install --prefix=/some/path

installation of C library
----------------------------

    git clone git@github.com:esheldon/psfex.git

    cd psfex

    # to install globally
    make install

    # to install in a particular place
    make install prefix=/some/path


dependencies
------------

- numpy
- cfitsio
