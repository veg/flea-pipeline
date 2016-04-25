#!/usr/bin/env python

"""
Get FASTA file with CCSs mapped to a single HQCS during the copynumber
task.

Usage:
  copynumber_to_fasta.py <hqcs_id> <pairs> <fasta> <outfile>
  copynumber_to_fasta.py -h | --help

Options:
  -h --help

"""

import os
from glob import glob
from itertools import chain
from collections import defaultdict

from Bio import SeqIO
from docopt import docopt


if __name__ == "__main__":
    args = docopt(__doc__)
    hqcs_id = args["<hqcs_id>"]
    pairfile = args["<pairs>"]
    fastafile = args["<fasta>"]
    outfile =  args["<outfile>"]

    with open(pairfile) as f:
        lines = list(line for line in f.read().split("\n") if line)
        pairs = list(line.split("\t") for line in lines)

    records = (SeqIO.parse(fastafile, 'fasta'))
    ccs_ids = set(a for a, b in pairs if b == hqcs_id)
    result = (r for r in records if r.id in ccs_ids)
    SeqIO.write(result, outfile, 'fasta')
