#!/usr/bin/env python

from setuptools import setup, find_packages
from flea_pipeline import __version__

setup(name='flea-pipeline',
      version=__version__,
      description='Full-length envelope analyzer',
      author='Kemal Eren',
      author_email='kemal@kemaleren.com',
      packages=find_packages(exclude="test"),
      test_suite='nose.collector',
      setup_requires=[
        'nose',
        'coverage',
        ],
      scripts=['flea_pipeline/flea.py'],
     )
