#!/usr/bin/env python
"""
Embedding based on precomputed distances.

"""
import json

import click
import numpy as np
import pandas as pd
from sklearn.manifold import MDS


def parse_file(infile):
    df = pd.read_csv(infile)
    labels = set()
    labels.update(df['ID1'])
    labels.update(df['ID2'])
    nseqs = len(labels)
    labels = list(labels)
    label_idx = dict((label, i) for i, label in enumerate(labels))
    X = np.zeros((nseqs, nseqs))
    for _, row in df.iterrows():
        X[label_idx[row['ID1']], label_idx[row['ID2']]] = row['Distance']
        X[label_idx[row['ID2']], label_idx[row['ID1']]] = row['Distance']
    return labels, X


@click.command()
@click.argument('infile')
@click.argument('outfile')
@click.option('--n-jobs', default=1)
def mds_cluster_file(infile, outfile, n_jobs):
    labels, X = parse_file(infile)
    model = MDS(dissimilarity='precomputed',
                n_init=16, max_iter=1000, n_jobs=n_jobs)
    coords = model.fit_transform(X)
    data = dict((label, list(c)) for label, c in zip(labels, coords))
    with open(outfile, 'w') as f:
        json.dump(data, f)


if __name__ == '__main__':
    mds_cluster_file()
