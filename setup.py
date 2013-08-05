import os
from distutils.core import setup, Extension
import numpy

ext=Extension("psfex._psfex_pywrap", 
              ["psfex/psfex_pywrap.c","psfex/psfex.c","psfex/poly.c"],
              extra_compile_args = ['-std=gnu99'])

setup(name="psfex", 
      version="0.1.0",
      description="Python and C libraries for reconstruct PSFEx psfs",
      license = "GPL",
      author="Erin Scott Sheldon",
      author_email="erin.sheldon@gmail.com",
      ext_modules=[ext],
      include_dirs=numpy.get_include(),
      packages=['psfex'])
