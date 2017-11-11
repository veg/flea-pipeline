#!/usr/bin/env python
"""
Embedding based on precomputed distances.

"""
import json

import click
import numpy as np
import pandas as pd
from sklearn.manifold import MDS

from flea.util import get_date_dict, id_to_label, id_to_copynumber


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


def weighted_median(values, weights):
    if len(values) != len(weights):
        raise Exception('lengths do not match')
    if not values:
        raise Exception('cannot take median of nothing')

    values = np.array(values)

    weights = np.array(weights)
    weights = weights / weights.sum()

    indices = np.argsort(values)

    values = values[indices]
    weights = weights[indices]

    weight_sum = 0
    result = values[0]
    for val, weight in zip(values, weights):
        weight_sum += weight
        if weight_sum < 0.5:
            result = val
        else:
            break
    return result


def weighted_median_both(values, weights):
    a = weighted_median(values, weights)
    b = weighted_median(list(reversed(values)), list(reversed(weights)))
    return np.mean([a, b])


@click.command()
@click.argument('infile')
@click.argument('metafile')
@click.argument('outfile')
@click.option('--n-jobs', default=1)
def mds_cluster_file(infile, metafile, outfile, n_jobs):
    ids, X = parse_file(infile)
    model = MDS(dissimilarity='precomputed',
                n_init=16, max_iter=1000, n_jobs=n_jobs)
    coords = model.fit_transform(X)

    # try to rotate to make time run left to right.
    #
    # could probably do something like PCA, but for now this works.
    #
    # finds the copynumber-weighted center of mass of the earliest and
    # latest time points, and finds the rotation that puts them both
    # on the x-axis, with the later time point to the right of the
    # earliest.
    date_dict = get_date_dict(metafile)
    dates = list(date_dict[id_to_label(id_)] for id_ in ids)
    cns = list(id_to_copynumber(id_) for id_ in ids)

    earliest_date = min(date_dict.values())
    latest_date = max(date_dict.values())

    early_xs = list(x for (x, y), date in zip(coords, dates) if date == earliest_date)
    early_ys = list(y for (x, y), date in zip(coords, dates) if date == earliest_date)
    early_cns = list(cn for cn, date in zip(cns, dates) if date == earliest_date)

    late_xs = list(x for (x, y), date in zip(coords, dates) if date == latest_date)
    late_ys = list(y for (x, y), date in zip(coords, dates) if date == latest_date)
    late_cns = list(cn for cn, date in zip(cns, dates) if date == latest_date)

    # use weighted median to find center of gravity, to downweight outliers
    early_x = weighted_median_both(early_xs, early_cns)
    early_y = weighted_median_both(early_ys, early_cns)
    late_x = weighted_median_both(late_xs, late_cns)
    late_y = weighted_median_both(late_ys, late_cns)

    x1 = late_x - early_x
    y1 = late_y - early_y
    x2 = 0
    y2 = -1

    rotation = np.array([[x1 * x2 + y1 * y2, -(x1 * y2 - x2 * y1)],
                         [x1 * y2 - x2 * y1, x1 * x2 + y1 * y2]])

    data = dict((id_, list(np.dot(rotation, c))) for id_, c in zip(ids, coords))
    with open(outfile, 'w') as f:
        json.dump(data, f)


if __name__ == '__main__':
    mds_cluster_file()
