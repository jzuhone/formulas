#!/usr/bin/env python
from setuptools import setup

setup(name='formulas',
      packages=['formulas'],
      version='0.1.0',
      description='Python tools for ACIS Ops',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/jzuhone/formulas',
      download_url='https://github.com/jzuhone/formulas/tarball/0.1.0',
      install_requires=["six","numpy","astropy","yt"],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      )
