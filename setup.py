# -*- encoding: UTF-8 -*-
from setuptools import setup, find_packages

with open('README.md','r') as f:
    long_description=f.read()

VERSION = '9.1_barcode_list'

setup(name='TECT',
      version=VERSION,
      description="Get the expression of TE from single cell RNA sequencing data, then Use the expression to improve clustering effect.",
      long_description=long_description,
      classifiers=[
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.8'
      ],
      keywords='python scRNA-seq TE cluster',
      author='yzfff6',
      author_email='yzfff6@outlook.com',
      url='https://github.com/yzf072/TECounter',
      packages=[
          'TECT'
      ],
      install_requires=[
        'sortedcontainers',
        'pysam',
        'numpy==1.19.5',
        'scanpy==1.7.0',
        'anndata==0.7.5',
        'umap-learn==0.5.1',
        'numpy==1.19.5',
        'scipy==1.6.0',
        'pandas==1.2.1',
        'scikit-learn==0.24.1'
      ],
      scripts=[
          'bin/TECT',
          'bin/TECT_clust',
          'bin/TECT_build'
      ]
      )
