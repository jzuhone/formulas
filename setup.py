#!/usr/bin/env python
from setuptools import setup
import sys
import os

if sys.argv[-1] == 'test':
    test_requirements = [
        'pytest',
        'coverage'
    ]
    try:
        modules = map(__import__, test_requirements)
    except ImportError as e:
        err_msg = e.message.replace("No module named ", "")
        msg = "%s is not installed. Install your test requirments." % err_msg
        raise ImportError(msg)
    os.system('py.test')
    sys.exit()

setup(name='formulas',
      packages=['formulas'],
      version='0.1.0',
      description='Unitful physical formulae',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/jzuhone/formulas',
      download_url='https://github.com/jzuhone/formulas/tarball/0.1.0',
      install_requires=["six","numpy","astropy","yt","matplotlib"],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      )