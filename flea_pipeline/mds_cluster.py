#!/usr/bin/env python
"""
MDS clustering based on precomputed distances.

Usage:
  mds_cluster.py [options] <infile> <outfile>
  mds_cluster.py -h | --help

Options:
  --flip        Flip values.
  -h --help     Show this screen.

"""
import json
import numpy as np
from docopt import docopt
from sklearn.manifold import MDS


def parse_file(infile):
    with open(infile) as f:
        lines = f.read().strip().split('\n')
    labels = set()
    for line in lines:
        label1, label2, val = line.split()
        labels.add(label1)
        labels.add(label2)
    nseqs = len(labels)
    labels = list(labels)
    label_idx = dict((label, i) for i, label in enumerate(labels))
    X = np.zeros((nseqs, nseqs))
    for line in lines:
        label1, label2, val = line.split()
        X[label_idx[label1], label_idx[label2]] = float(val)
        X[label_idx[label2], label_idx[label1]] = float(val)
    return labels, X


def mds_cluster_file(infile, outfile, do_flip):
    labels, X = parse_file(infile)
    if do_flip:
        # max value should correspond to totally similar items
        X = X.max() - X
    model = MDS(dissimilarity='precomputed')
    coords = model.fit_transform(X)
    data = dict((label, list(c)) for label, c in zip(labels, coords))
    with open(outfile, 'w') as f:
        json.dump(data, f)


if __name__ == '__main__':
    args = docopt(__doc__)
    infile = args["<infile>"]
    outfile = args["<outfile>"]
    do_flip = args['--flip']
    mds_cluster_file(infile, outfile, do_flip)
