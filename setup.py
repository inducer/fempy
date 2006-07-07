#!/usr/bin/env python
# -*- coding: latin-1 -*-

from distutils.core import setup,Extension
import glob
import os
import os.path

setup(name="fempy",
      version="0.10",
      description="A finite element package",
      author=u"Andreas Klöckner",
      author_email="inform@tiker.net",
      license = "GNU GPL",
      url="http://news.tiker.net/software/fempy",
      packages=["fempy"],
      package_dir={"fempy": "src"}
     )
