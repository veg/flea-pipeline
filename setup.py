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
        'flea_pipeline/trim_terminal_gaps.py',
        'flea_pipeline/diagnose.py',
        'flea_pipeline/DNAcons.py',
        'flea_pipeline/perfectORFs.py',
        'flea_pipeline/translate.py',
        'flea_pipeline/cluster_fastq.py',
        'flea_pipeline/trim.py',
        'flea_pipeline/mds_cluster.py',
        'scripts/publish.py',
        'scripts/parse_log.py',
        'flea_pipeline/HXB2partsSplitter.bf',
        'flea_pipeline/obtainEvolutionaryHistory.bf',
        'flea_pipeline/runFUBAR.bf',
        'flea_pipeline/dates.bf',
        'flea_pipeline/tools.bf',
        ],
      )
