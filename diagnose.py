#!/usr/bin/env python

"""
Compare HQCS frequencies to CCS frequencies.

Usage:
  diagnose.py [options] <hqcs_file> <ccs_file> <copynumber_file> <result_path>
  diagnose.py -h | --help

Options:
  -c <FLOAT>  Cutoff for JS distance [default: 0.01]
  -h --help   Show this screen

"""


import csv
import os
import warnings

from docopt import docopt
import numpy as np
from Bio import SeqIO

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)


def to_freqs(a):
    return a / a.sum(axis=0)


def kl_divergence(p, q):
    if np.any(q == 0):
        raise Exception('cannot compute KL divergence because q == 0')
    with warnings.catch_warnings():
        # can safely ignore warnings here; we will filter out invalid values
        warnings.simplefilter("ignore")
        elts = (p * (np.log(p) - np.log(q)))
    elts[p == 0] = 0
    result = elts.sum()
    if np.isnan(result):
        raise Exception('KL divergence failed')
    return result


def js_divergence(p, q):
    # ensure there are no zeros
    m = (p + q) / 2
    bools = (m > 0)
    p = p[bools]
    q = q[bools]
    m = m[bools]
    return (kl_divergence(p, m) + kl_divergence(q, m)) / 2


def column_count(a, keys, weights=None):
    """
    >>> column_count(np.array([['A', 'A'], ['A', 'B']]), ['A', 'B'])
    array([[2, 1],
           [0, 1]])

    >>> column_count(np.array([['A', 'A'], ['A', 'B']]), keys=['A', 'B'], weights=[3, 2])
    array([[5, 3],
           [0, 2]])

    """
    keys = np.array(keys).ravel()
    result = (a == keys[:, np.newaxis, np.newaxis])
    if weights is not None:
        weights = np.array(weights).ravel()
        assert len(weights) == len(a)
        result = result * weights.reshape(-1, 1)
        assert np.all(result.sum(axis=1).sum(axis=0) == weights.sum())
    return result.sum(axis=1)


def diagnose(filename, filename_ccs, copynumbers_file, result_path, cutoff):
    def out(f):
        return os.path.join(result_path, f)

    format = "fasta"

    #Getting and sorting copynumbers
    with open(copynumbers_file, 'r') as handle:
        parsed = csv.reader(handle, delimiter='\t')
        cns = list(sorted(parsed))
        copynumber_array = np.array(list(int(n) for elem, n in cns))

    #Getting and sorting HQCS seqs
    hqcs = list(sorted(SeqIO.parse(filename, format), key=lambda r: r.id))
    hqcs_array = np.array(list(list(str(rec.seq)) for rec in hqcs))
    hqcs_labels = np.array(list(rec.id.split("_")[0] for rec in hqcs))

    if len(copynumber_array) != len(hqcs_array):
        raise Exception('number of HQCS sequence does not match copy numbers')

    # import ccs seqs
    ccs = list(SeqIO.parse(filename_ccs, format))
    ccs_array = np.array(list(list(str(rec.seq)) for rec in ccs))
    ccs_labels = np.array(list(rec.id.split("_")[0] for rec in ccs))

    if hqcs_array.shape[1] != ccs_array.shape[1]:
        raise Exception('alignment shapes do not match')

    # amino acids
    aas = list(sorted(set(hqcs_array.ravel()) | set(ccs_array.ravel())))
    keys = np.array(aas)

    # Relying on "_" to separate timepoint ID
    labels = sorted(set(hqcs_labels))

    for label in labels:
        bools = hqcs_labels == label
        hqcs_counts = column_count(hqcs_array[bools], keys, weights=copynumber_array[bools])

        ccs_counts = column_count(ccs_array[ccs_labels == label], keys)

        np.savetxt(out("hqcs_counts_{}.csv".format(label)), hqcs_counts.T,
                   fmt="%1.1u", delimiter=",", header=",".join(aas))
        np.savetxt(out("ccs_counts_{}.csv".format(label)), ccs_counts.T,
                   fmt="%u", delimiter=",", header=",".join(aas))

        # plot all freqs
        plt.scatter(np.ravel(to_freqs(hqcs_counts)),
                    np.ravel(to_freqs(ccs_counts)),
                    alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('CCS frequency')
        plt.title(label)
        plt.savefig(out("freq_agreement_{}.png".format(label)))
        plt.close()

        # strip X and normalize
        x_index = aas.index('X')
        hqcs_no_x = np.delete(hqcs_counts, x_index, 0)
        hqcs_no_x_freqs = to_freqs(hqcs_no_x)
        ccs_no_x = np.delete(ccs_counts, x_index, 0)
        ccs_no_x_freqs = to_freqs(ccs_no_x)

        # plot without X
        plt.scatter(np.ravel(hqcs_no_x_freqs), np.ravel(ccs_no_x_freqs), alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('CCS frequency')
        plt.title("{} no X".format(label))
        plt.savefig(out("freq_agreement_no_x_{}.png".format(label)))
        plt.close()

        divs = np.array(list(js_divergence(a, b)
                             for a, b in zip(hqcs_no_x_freqs.T, ccs_no_x_freqs.T)))

        bad_columns = np.where(divs > cutoff)[0]

        # plot column-wise JS divergence
        plt.plot(divs)
        if len(bad_columns) > 0:
            plt.axhline(y=cutoff, color='r')
        plt.xlabel('column')
        plt.ylabel('JS divergence')
        plt.title("{} no X".format(label))
        plt.savefig(out("js_div_no_x_{}.png".format(label)))
        plt.close()

        # compute bad locations which have extreme JS divergence

        if len(bad_columns) > 0:
            cols = [bad_columns,
                    divs[bad_columns],
                    keys[hqcs_no_x_freqs[:, bad_columns].argmax(axis=0)],
                    hqcs_no_x_freqs[:, bad_columns].max(axis=0),
                    keys[ccs_no_x_freqs[:, bad_columns].argmax(axis=0)],
                    ccs_no_x_freqs[:, bad_columns].max(axis=0)]
            rows = list(zip(*cols))
            with open(out("bad_sites_{}.csv".format(label)), "w") as f:
                writer = csv.writer(f)
                writer.writerow(['column', 'js_divergence', 'hqcs_aa', 'hqcs_freq', 'ccs_aa', 'ccs_freq'])
                writer.writerows(rows)


if __name__ == "__main__":
    args = docopt(__doc__)
    hqcs_file = args["<hqcs_file>"]
    ccs_file = args["<ccs_file>"]
    copy_file = args["<copynumber_file>"]
    result_path = args["<result_path>"]
    cutoff = float(args["-c"])
    diagnose(hqcs_file, ccs_file, copy_file, result_path, cutoff)
