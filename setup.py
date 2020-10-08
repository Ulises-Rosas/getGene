#!/usr/bin/env python3

import setuptools
from distutils.core import setup

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(name="genetable",
      version='0.1',
      author='Ulises Rosas',
      long_description = readme,
      long_description_content_type = 'text/markdown',
      author_email='ulisesfrosasp@gmail.com',
      url='https://github.com/Ulises-Rosas/geneTable',
      packages = ['genetable'],
      package_dir = {'genetable': 'src'},
      entry_points={
        'console_scripts': [
            'getgenomes  = genetable.getGenomes:main',
            'getFeatures = genetable.geneTable:main',
            'looksra     = genetable.lookSRA:main'
            ]
      },
      install_requires=['matplotlib'],
      classifiers = [
             'Programming Language :: Python :: 3',
             'License :: OSI Approved :: MIT License'
             ]
    )
