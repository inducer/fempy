#!/usr/bin/env python

from distutils.core import setup,Extension
import glob
import os
import os.path

setup(name="femPy",
      version="0.10",
      description="A finite element package",
      author="Andreas Kloeckner",
      author_email="ak@ixion.net",
      license = "GNU GPL",
      url="http://fempy.sf.net",
      packages=["fempy"],
      package_dir={"fempy": "src"}
     )
