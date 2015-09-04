#!/usr/bin/env python

"""
Compare HQCS frequencies to CCS frequencies.

Usage:
  diagnose.py [options] <hqcs_file> <ccs_file> <copynumber_file> <result_path>
  diagnose.py -h | --help

Options:
  -h --help               Show this screen

"""


import csv
import os

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
    return (p * (np.log(p) - np.log(q))).sum()


def js_divergence(p, q):
    m = 0.5 * (p + q)
    return 0.5 * (kl_divergence(p, m) + kl_divergence(q, m))


def column_count(a, keys, weights=None):
    keys = np.array(keys).ravel()
    result = (a == keys[:, np.newaxis, np.newaxis])
    if weights is not None:
        weights = weights.ravel()
        assert len(weights) == len(a)
        result = result * weights.reshape(-1, 1)
        assert np.all(result.sum(axis=1).sum(axis=0) == weights.sum())
    return result.sum(axis=1)


def diagnose(filename, filename_ccs, copynumbers_file, result_path):
    def out(f):
        return os.path.join(result_path, f)

    format = "fasta"

    #Getting and sorting copynumbers
    with open(copynumbers_file, 'r') as handle:
        parsed = csv.reader(handle, delimiter='\t')
        copynumber_array = np.array(list(int(n) for elem, n in sorted(parsed)))

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

    # TODO: get a general threshold region for bad pairs, and use this instead.
    thresh = 3
    epsil = 0.02

    for label in labels:
        bools = hqcs_labels == label
        hqcs_counts = column_count(hqcs_array[bools], keys, weights=copynumber_array[bools])
        hqcs_freqs = to_freqs(hqcs_counts)

        ccs_counts = column_count(ccs_array[ccs_labels == label], keys)
        ccs_freqs = to_freqs(ccs_counts)

        np.savetxt(out("hqcs_counts_{}.csv".format(label)), hqcs_counts.T,
                   fmt="%1.1u", delimiter=",", header=",".join(aas))
        np.savetxt(out("ccs_counts_{}.csv".format(label)), ccs_counts.T,
                   fmt="%u", delimiter=",", header=",".join(aas))

        # compute bad locations, which have extreme HQCS counts, but less
        # extreme CCS counts.

        # Don't ask about this parameterization. It separated out the ones I needed to see.
        bad_inds = np.greater(ccs_freqs, -epsil + epsil * thresh + thresh * hqcs_freqs)
        bad_inds_rev = np.greater(1 - ccs_freqs, -epsil + epsil * thresh + thresh * (1 - hqcs_freqs))
        combined = np.logical_or(bad_inds, bad_inds_rev)

        bad_keys, bad_columns = np.where(combined)
        rows = zip(*[bad_columns, keys[bad_keys],
                     hqcs_freqs[bad_keys, bad_columns],
                     ccs_freqs[bad_keys, bad_columns]])
        with open(out("bad_sites_{}.csv".format(label)), "w") as f:
            writer = csv.writer(f)
            writer.writerow(['column', 'aa', 'hqcs_freq', 'ccs_freq'])
            writer.writerows(rows)

        # plot all freqs
        plt.scatter(np.ravel(hqcs_freqs), np.ravel(ccs_freqs), alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('CCS frequency')
        plt.title(label)
        plt.savefig(out("freq_agreement_{}.png".format(label)))
        plt.close()

        # strip X and normalize
        x_index = aas.index('X')
        hqcs_no_x = np.delete(hqcs_counts, x_index, 0)
        ccs_no_x = np.delete(ccs_counts, x_index, 0)

        # plot without X
        plt.scatter(np.ravel(to_freqs(hqcs_no_x)), np.ravel(to_freqs(ccs_no_x)), alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('CCS frequency')
        plt.title("{} no X".format(label))
        plt.savefig(out("freq_agreement_no_x_{}.png".format(label)))
        plt.close()

        divs = list(js_divergence((a + 1) / (a.sum() + 1), (b + 1) / (b.sum() + 1))
                    for a, b in zip(hqcs_no_x.T, ccs_no_x.T))

        # plot column-wise JS divergence
        plt.plot(divs)
        plt.xlabel('column')
        plt.ylabel('JS divergence')
        plt.title("{} no X".format(label))
        plt.savefig(out("js_div_no_x_{}.png".format(label)))
        plt.close()


if __name__ == "__main__":
    args = docopt(__doc__)
    hqcs_file = args["<hqcs_file>"]
    ccs_file = args["<ccs_file>"]
    copy_file = args["<copynumber_file>"]
    result_path = args["<result_path>"]
    diagnose(hqcs_file, ccs_file, copy_file, result_path)
