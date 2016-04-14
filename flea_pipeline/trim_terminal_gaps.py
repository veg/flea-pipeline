#!/usr/bin/env python

"""Reverse complement (if necessary) and trim terminal gaps in a fastq
file.

<userfile> should have the following tab-seperated information, one
line per sequence, as in the usearch -userout file:
qstrand+tstrand+caln

Usage:
  trim_terminal_gaps.py [options] <fqfile> <userfile> <outfile>
  trim_terminal_gaps.py -h

Options:
  -v --verbose      Print summary [default: False]
  -h --help         Show this screen

"""

import csv
import re

from docopt import docopt

from Bio import SeqIO


pattern = r'[0-9]*[IDM]'


def n_del(part):
    """

    >>> n_del("D")
    1

    >>> n_del("34D")
    34

    >>> n_del("I")
    0

    >>> n_del("3I")
    0

    """
    if not re.match(pattern, part):
        raise Exception('argument does not match regexp')
    if len(part) == 1:
        return 1 if part == 'D' else 0
    return int(part[:-1]) if part[-1] == 'D' else 0


def n_trim(caln):
    """Get amount to trim from start and end.

    >>> n_trim("3D4M2D")
    (3, 2)

    >>> n_trim("3I4M2D")
    (0, 2)

    >>> n_trim("3I4M")
    (0, 0)

    """
    parts = re.findall(pattern, caln)
    return n_del(parts[0]), n_del(parts[-1])


def handle_record(record, query_rcomp, db_rcomp, caln):
    if query_rcomp != db_rcomp:
        rcomp = record.reverse_complement()
        rcomp.id = record.id
        rcomp.description = record.description
        record = rcomp
    a, b = n_trim(caln)
    start = a
    stop = len(record) - b
    return record[start:stop]


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<fqfile>"]
    userfile = args["<userfile>"]
    outfile = args["<outfile>"]

    records = SeqIO.parse(infile, 'fastq')
    with open(userfile) as handle:
        lines = list(line.split('\t') for line in handle.read().split('\n'))

    result = list(handle_record(r, qstrand, tstrand, caln)
                  for r, (qstrand, tstrand, caln) in zip(records, lines))

    SeqIO.write(result, outfile, 'fastq')
