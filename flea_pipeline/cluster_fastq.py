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
    """Returns a dict from cluster number to labels.

    Assumes input was generated with -top_hit_only.

    """
    result = defaultdict(list)
    with open(infile) as h:
        lines = h.read().strip().split('\n')
    records = list(line.split('\t') for line in lines)
    for elts in records:
        if elts[0] in 'SH':
            cluster = elts[1]
            label = elts[8]
            result[cluster].append(label)
    for k, v in result.items():
        if len(v) != len(set(v)):
            raise Exception('cluster {} contains duplicate ids'.format(k))
    for elts in records:
        if elts[0] == 'C':
            cluster = elts[1]
            size = int(elts[2])
            obs = len(result[cluster])
            if obs != size:
                raise Exception('cluster {} size {} != {}'.format(cluster, size, obs))
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
        outfile = os.path.join(outdir, "cluster_{}.raw.fastq".format(cluster_id))
        SeqIO.write(cluster_records, outfile, 'fastq')
