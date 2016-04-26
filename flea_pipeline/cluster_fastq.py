#!/usr/bin/env python

"""Write each cluster from a .uc file to a fastq file.

Each output file is put in <outdir> and named "cluster_{}.fastq".

Usage:
  cluster_fastq.py [options] <ucfile> <fastqfile> <outdir>
  cluster_fastq.py -h

Options:
  --minsize <INT>   Minimum cluster size [default: 1]
  -v --verbose      Print summary [default: False]
  -h --help         Show this screen

"""
import os
from collections import defaultdict

from docopt import docopt

from Bio import SeqIO


def parse_ucfile(infile):
    """Returns a dict of from cluster number to labels"""
    result = defaultdict(list)
    with open(infile) as h:
        lines = h.read().split('\n')
    for elts in (line.split('\t') for line in lines if line):
        cluster_id = elts[1]
        label = elts[8]
        result[cluster_id].append(label)
    return result


if __name__ == "__main__":
    args = docopt(__doc__)
    ucfile = args["<ucfile>"]
    fastqfile = args["<fastqfile>"]
    outdir = args["<outdir>"]
    minsize = int(args["--minsize"])

    records = SeqIO.parse(fastqfile, 'fastq')
    rdict = dict((r.id, r) for r in records)
    cdict = parse_ucfile(ucfile)

    for cluster_id, labels in cdict.items():
        cluster_records = list(rdict[label] for label in labels)
        if len(cluster_records) < minsize:
            continue
        outfile = os.path.join(outdir, "cluster_{}.fastq".format(cluster_id))
        SeqIO.write(cluster_records, outfile, 'fastq')
