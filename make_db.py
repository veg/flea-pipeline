#!/usr/bin/env python
"""
Makes a usearch database.

Usage:
  make_db.py <infile> <dbfile>
  make_db.py -h | --help

Options:
  -h --help      Show this screen.

"""
import sys
import subprocess
import shlex

from docopt import docopt

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped


def make_db(infile, dbfile):
    """Make usearch database of reads in `infile`."""
    kwargs = dict(infile=infile, dbfile=dbfile)
    cmd = "usearch -makeudb_usearch {infile} -output {dbfile}".format(**kwargs)
    subprocess.check_call(shlex.split(cmd))


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    dbfile = args["<dbfile>"]
    make_db(infile, dbfile)
