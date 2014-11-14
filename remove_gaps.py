#!/usr/bin/env python
"""
Removes gaps from MSA fasta to prepare clean reads for
align_to_clean.py.

Usage:
  remove_gaps.py <clean_msa_fasta>
  remove_gaps.py -h | --help

Options:
  -h --help      Show this screen.

"""
import sys

from docopt import docopt
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped


def ungap(record):
    record.seq = record.seq.ungap()
    return record


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<clean_msa_fasta>"]
    handle = open(infile, "rU")
    records = SeqIO.parse(handle, "fasta", alphabet=Gapped(IUPAC.unambiguous_dna))
    results = (ungap(r) for r in records)
    SeqIO.write(results, sys.stdout, "fasta")
