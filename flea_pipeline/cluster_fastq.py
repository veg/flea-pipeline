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
import random
from collections import defaultdict

from docopt import docopt

from Bio import SeqIO


def parse_ucfile(infile):
    """Returns a dict from cluster number to labels.

    Assumes input was generated with -top_hit_only. Breaks tie at
    random, so cluster sizes may not exactly match those reported by
    usearch.

    """
    centroid_to_cluster = {}
    seq_to_centroids = defaultdict(list)
    with open(infile) as h:
        lines = h.read().strip().split('\n')
    records = list(line.split('\t') for line in lines)
    for elts in records:
        if elts[0] == 'S':
            cluster = elts[1]
            centroid = elts[8]
            centroid_to_cluster[centroid] = cluster
        elif elts[0] == 'H':
            label = elts[8]
            centroid = elts[9]
            seq_to_centroids[label].append(centroid)
    result = defaultdict(list)
    for label, centroids in seq_to_centroids.items():
        centroid = random.choice(centroids)
        result[centroid_to_cluster[centroid]].append(label)

    for k, v in result.items():
        if len(v) != len(set(v)):
            raise Exception('cluster {} contains duplicate ids'.format(k))
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
