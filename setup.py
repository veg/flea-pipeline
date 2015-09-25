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
      scripts=[
        'flea_pipeline/flea.py',
        'flea_pipeline/backtranslate.py',
        'flea_pipeline/correct_shifts.py',
        'flea_pipeline/diagnose.py',
        'flea_pipeline/DNAcons.py',
        'flea_pipeline/perfectORFs.py',
        'flea_pipeline/translate.py',
        'flea_pipeline/trim.py',
        'scripts/publish.py',
        ],
     )
