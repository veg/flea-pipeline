#!/usr/bin/env python

from setuptools import setup, find_packages
from flea import __version__

setup(name='flea',
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
        'flea/backtranslate.py',
        'flea/correct_shifts.py',
        'flea/trim_terminal_gaps.py',
        'flea/diagnose.py',
        'flea/DNAcons.py',
        'flea/perfectORFs.py',
        'flea/translate.py',
        'flea/cluster_fastq.py',
        'flea/trim.py',
        'flea/manifold_embed.py',
        'scripts/publish.py',
        ],
      )
