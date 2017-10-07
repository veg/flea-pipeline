#!/usr/bin/env python

from setuptools import setup, find_packages
from flea import __version__

setup(name='flea',
      version=__version__,
      description='Helper scripts for FLEA pipeline.',
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
        'flea/cluster_fastq.py',
        'flea/coordinates_json.py',
        'flea/correct_shifts.py',
        'flea/diagnose.py',
        'flea/DNAcons.py',
        'flea/filter_fastx.py',
        'flea/insert_gaps.py',
        'flea/js_divergence.py',
        'flea/manifold_embed.py',
        'flea/pairs_to_fasta.py',
        'flea/sequences_json.py',
        'flea/translate.py',
        'flea/trim_tails.py',
        'flea/write_copynumbers.py',
        'scripts/publish.py',
        ],
      )
