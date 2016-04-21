#!/usr/bin/env python

"""
Map each HQCS id to all CCS ids that generated it.

There are multiple CCSs per cluster, and multiple clusters per HQCS,
since unique shift-corrected HQCSs get combined.

Takes as input an "alignment" directory from a pipeline run.

Usage:
  hqcs_to_ccs_ids.py [options] <directory>
  hqcs_to_ccs_ids.py | --help

Options:
  -o <outdir>  Output directory
  -h --help

"""

import os
from glob import glob
from itertools import chain
from collections import defaultdict

from Bio import SeqIO
from docopt import docopt


def hqcs_to_ccs(directory):
    # for each *.clusters/*raw.fasta, get CCS ids
    # map to id in *.aligned.hqcs.fasta
    files = glob(os.path.join(directory, "*.clusters/*.raw.fasta"))
    result = {}
    for f in files:
        ids = list(r.id for r in SeqIO.parse(f, 'fasta'))
        hqcs_file = ''.join([os.path.splitext(f)[0], '.keep.aligned.hqcs.fasta'])
        hqcs = list(SeqIO.parse(hqcs_file, 'fasta'))
        if hqcs:
            id_ = hqcs[0].id
            result[id_] = ids
    return result


def map_unique(src_recs, unique_recs):
    """result[<unique_id>] = list of src ids with identical sequences"""
    result = defaultdict(list)
    seq_to_id = dict((str(r.seq), r.id) for r in unique_recs)
    for r in src_recs:
        id_ = seq_to_id[str(r.seq)]
        result[id_].append(r.id)
    return result


def unique_to_hqcs(directory):
    result = {}
    src_files = glob(os.path.join(directory, "*.shifted.fasta"))
    for f in src_files:
        unique = ''.join([os.path.splitext(f)[0], '.uniques.fasta'])
        src_recs = list(SeqIO.parse(f, 'fasta'))
        unique_recs = list(SeqIO.parse(unique, 'fasta'))
        d = map_unique(src_recs, unique_recs)
        for k in d:
            assert k not in result
        result.update(d)
    return d


def unique_to_ccs(hqcs_to_ccs_ids, unique_to_hqcs_ids):
    """map every unique HQCS's id to all the CCS ids in the clusters
    that generated it.

    """
    result = {}
    for id_, ids in unique_to_hqcs_ids.items():
        ccs_ids = list(hqcs_to_ccs_ids[i] for i in ids)
        result[id_] = list(chain(*ccs_ids))
    return result


if __name__ == "__main__":
    args = docopt(__doc__)
    directory = args["<directory>"]
    if args["-o"]:
        outdir = args["-o"]
    else:
        outdir = directory

    hqcs_to_ccs_ids = hqcs_to_ccs(directory)
    unique_to_hqcs_ids = unique_to_hqcs(directory)

    unique_to_ccs_ids = unique_to_ccs(hqcs_to_ccs_ids, unique_to_hqcs_ids)

    outfile = os.path.join(outdir, "unique-to-ccs-ids.tsv")
    with open(outfile, 'w') as h:
        for key, value in unique_to_ccs_ids.items():
            h.write("{}\t{}\n".format(key, '\t'.join(value)))

    outfile = os.path.join(outdir, "unique-to-ccs-counts.tsv")
    with open(outfile, 'w') as h:
        for key, value in unique_to_ccs_ids.items():
            h.write("{}\t{}\n".format(key, len(value)))

    outfile = os.path.join(outdir, "unique-to-hqcs-ids.tsv")
    with open(outfile, 'w') as h:
        for key, value in unique_to_hqcs_ids.items():
            h.write("{}\t{}\n".format(key, '\t'.join(value)))

    outfile = os.path.join(outdir, "unique-to-hqcs-counts.tsv")
    with open(outfile, 'w') as h:
        for key, value in unique_to_hqcs_ids.items():
            h.write("{}\t{}\n".format(key, len(value)))
