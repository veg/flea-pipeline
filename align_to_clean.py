#!/usr/bin/env python
"""
Aligns each read to its best-matched clean read.

Dependencies:
  - usearch
  - balign

Usage:
  align_to_clean.py <clean_fasta> <dirty_fasta>
  align_to_clean.py -h | --help

Options:
  --db <dbfile>  Where to build usearch db. Defaults to /tmp/{rand}.udb
  --keep         Do not delete dbfile after completion.
  -h --help      Show this screen.

"""
import itertools
import subprocess
import shlex
import uuid
import os

from docopt import docopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


# Taken from itertools recipes
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks."
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def make_db(infile, dbfile):
    kwargs = dict(infile=infile, dbfile=dbfile)
    cmd = "usearch -makeudb_usearch {infile} -output {dbfile}".format(**kwargs)
    subprocess.check_call(shlex.split(cmd))


def usearch_global(readfile, dbfile, pairfile):
    kwargs = dict(readfile=readfile, dbfile=dbfile, pairfile=pairfile)
    cmd = ("usearch -usearch_global {readfile} -db {dbfile} -id 0.8"
           "-userout {pairfile} -userfields query+target -strand both")
    cmd = cmd.format(**kwargs)
    subprocess.check_call(shlex.split(cmd))


def match_reads(clean_filename, dirty_filename):
    """For each dirty read, find the best clean read.

    Uses usearch to build a database of clean reads and run global
    alignment.
    
    Returns: list of (dirty id, clean id) tuples.

    """
    tmpname = uuid.uuid4()
    dbfile = "/tmp/{}.udb".format(tmpname)
    pairfile = "/tmp/{}.fasta".format(tmpname)
    make_db(clean_filename, dbfile)
    usearch_global(dirty_filename, dbfile, pairfile)
    with open(pairfile) as f:
        result = list(line.split("\t") for line in f.readlines())
    os.rm(dbfile)
    os.rm(pairfile)
    return result


if __name__ == "__main__":
    args = docopt(__doc__)
    clean_filename = args["<clean_fasta>"]
    dirty_filename = args["<dirty_fasta>"]
    result = match_reads(clean_filename, dirty_filename)
