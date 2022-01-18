#!/usr/bin/env python3

# import sys
# import platform
import setuptools
from distutils.core import setup

dependencies = ['requests']

# myos = sys.platform

# if myos == 'darwin':
#     bins = ['./ext_bin/darwin/vdb-config',
#             './ext_bin/darwin/prefetch',
#             './ext_bin/darwin/fasterq-dump']

# else:
#     sys.stderr.write('Currently genetable does not work with %s operative system'  % myos)
#     sys.stderr.flush()
#     exit()

with open('README.md') as readme_file:
    readme = readme_file.read()

# ext_bin_data = [
#     'darwin_fasterq-dump',
#     'darwin_prefetch',
#     'darwin_vdb-config',
#     'test.txt'
#     ]

setup(name="genetable",
      version='0.6',
      author='Ulises Rosas',
      long_description = readme,
      long_description_content_type = 'text/markdown',
      author_email='ulisesfrosasp@gmail.com',
      url='https://github.com/Ulises-Rosas/geneTable',
      packages = ['genetable'],
      package_dir = {'genetable': 'src'},
      entry_points={
        'console_scripts': [
            'lookfeatures = genetable.geneTable:main',
            'lookgenomes = genetable.lookgenomes:main',
            'looksra     = genetable.lookSRA:main',
            'getgenomes  = genetable.getgenomes:main',
            ]
      },
      install_requires=dependencies,
      classifiers = [
             'Programming Language :: Python :: 3',
             'License :: OSI Approved :: MIT License'
             ]
    )
