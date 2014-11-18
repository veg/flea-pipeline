#!/usr/bin/env python
"""
Removes gaps from reads and makes a usearch database.

Usage:
  make_db.py <infile> <outfile> <dbfile>
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


def ungap(record):
    record.seq = record.seq.ungap()
    return record


def remove_gaps(infile, outfile):
    records = SeqIO.parse(infile, "fasta", alphabet=Gapped(IUPAC.unambiguous_dna))
    results = (ungap(r) for r in records)
    # need to have length multiple of 3, for codon alignment downstream
    # TODO: warn if any have wrong length.
    results = (r for r in results if (len(r) % 3) == 0)
    SeqIO.write(results, outfile, "fasta")


def make_db(infile, dbfile):
    """Make usearch database of reads in `infile`."""
    kwargs = dict(infile=infile, dbfile=dbfile)
    cmd = "usearch -makeudb_usearch {infile} -output {dbfile}".format(**kwargs)
    subprocess.check_call(shlex.split(cmd))


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    outfile = args["<outfile>"]
    dbfile = args["<dbfile>"]
    remove_gaps(infile, outfile)
    make_db(outfile, dbfile)
