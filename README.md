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

    # if you install in a prefix, make sure you
    # add the /some/path/lib/python2.7/site-packages
    # directory to your PYTHONPATH

installation of C library
----------------------------

    git clone git@github.com:esheldon/psfex.git

    cd psfex

    # to install globally
    make install

    # to install in a particular place
    make install prefix=/some/path

    # if you install in a prefix, make sure you
    # add the /some/path/lib
    # directory to your LD_LIBRARY_PATH and
    # LIBRARY_PATH.  Similarly add /some/path/include
    # to C_LIBRARY_PATH and CPATH


dependencies
------------

- numpy
- cfitsio
