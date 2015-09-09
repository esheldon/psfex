import os
from distutils.core import setup, Extension
import numpy

ext=Extension("psfex._psfex_pywrap", 
              ["psfex/psfex_pywrap.c","psfex/psfex.c","psfex/poly.c"],
              extra_compile_args=['-std=gnu99'],
              include_dirs=[numpy.get_include()])

exec(open('psfex/version.py').read())

setup(name="psfex", 
      version=__version__,
      description="Python and C libraries for reconstruct PSFEx psfs",
      license = "GPL",
      author="Erin Scott Sheldon",
      author_email="erin.sheldon@gmail.com",
      ext_modules=[ext],
      packages=['psfex'])
