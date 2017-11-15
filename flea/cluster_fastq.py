#!/usr/bin/env python

"""Write each cluster from a .uc file to a FASTQ or FASTA file.

Each output file is put in <outdir> and named "cluster_{}_raw.fastq".

"""
import sys
import os
from collections import defaultdict

import click

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


@click.command()
@click.argument('ucfile')
@click.argument('outdir')
@click.option('--minsize', type=int, default=1, help="minimum cluster size")
@click.option('--fasta', is_flag=True, help='write to FASTA instead of FASTQ')
@click.option('-v', '--verbose', is_flag=True)
def cluster_fastq(ucfile, outdir, minsize, fasta, verbose):
    records = list(SeqIO.parse(sys.stdin, 'fastq'))
    rdict = dict((r.id, r) for r in records)
    cdict = parse_ucfile(ucfile)

    format = 'fasta' if fasta else 'fastq'

    for cluster_id, labels in cdict.items():
        cluster_records = list(rdict[label] for label in labels)
        if len(cluster_records) < minsize:
            continue
        outfile = os.path.join(outdir, "cluster_{}_raw.{}".format(cluster_id, format))
        SeqIO.write(cluster_records, outfile, format)


if __name__ == "__main__":
    cluster_fastq()
