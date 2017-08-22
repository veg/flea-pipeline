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

from docopt import docopt
import numpy as np
from Bio import SeqIO

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from flea_pipeline.util import column_count, prob, js_divergence


np.set_printoptions(suppress=True)


def diagnose(filename, filename_ccs, copynumbers_file, result_path, cutoff):
    def out(f):
        return os.path.join(result_path, f)

    format = "fasta"

    plt.style.use('seaborn-deep')

    #Getting and sorting copynumbers
    with open(copynumbers_file, 'r') as handle:
        parsed = csv.reader(handle, delimiter='\t')
        cns = list(sorted(parsed))
        copynumber_array = np.array(list(int(n) for _, n in cns))
        copynumber_ids = np.array(list(elem for elem, _ in cns))

    #Getting and sorting HQCS seqs
    hqcs = list(sorted(SeqIO.parse(filename, format), key=lambda r: r.id))
    hqcs_array = np.array(list(list(str(rec.seq)) for rec in hqcs))
    hqcs_labels = np.array(list(rec.id.split("_")[0] for rec in hqcs))
    hqcs_ids = np.array(list(rec.id for rec in hqcs))

    if len(copynumber_array) != len(hqcs_array):
        raise Exception('number of HQCS sequence does not match copy numbers')

    if not np.all(copynumber_ids == hqcs_ids):
        raise Exception('HQCS ids do not match copy number ids')

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

    for label in labels + ['overall']:
        if label == 'overall':
            hqcs_bools = np.array([True] * len(hqcs_array))
            ccs_bools = np.array([True] * len(ccs_array))
        else:
            hqcs_bools = hqcs_labels == label
            ccs_bools = ccs_labels == label

        hqcs_counts = column_count(hqcs_array[hqcs_bools], keys, weights=copynumber_array[hqcs_bools])
        ccs_counts = column_count(ccs_array[ccs_bools], keys)

        if np.all(hqcs_counts.sum(axis=0) != copynumber_array[hqcs_bools].sum()):
            raise Exception('bad hqcs counts')

        if np.all(ccs_counts.sum(axis=0) != ccs_bools.sum()):
            raise Exception('bad ccs counts')

        np.savetxt(out("hqcs_counts_{}.csv".format(label)), hqcs_counts.T,
                   fmt="%1.1u", delimiter=",", header=",".join(aas))
        np.savetxt(out("ccs_counts_{}.csv".format(label)), ccs_counts.T,
                   fmt="%u", delimiter=",", header=",".join(aas))

        # plot all freqs
        plt.scatter(np.ravel(prob(hqcs_counts, axis=0)),
                    np.ravel(prob(ccs_counts, axis=0)),
                    alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('CCS frequency')
        plt.title(label)
        plt.savefig(out("freq_agreement_{}.pdf".format(label)))
        plt.close()

        # strip X and normalize
        x_index = aas.index('X')
        hqcs_no_x = np.delete(hqcs_counts, x_index, 0)
        hqcs_no_x_freqs = prob(hqcs_no_x, axis=0)
        ccs_no_x = np.delete(ccs_counts, x_index, 0)
        ccs_no_x_freqs = prob(ccs_no_x, axis=0)

        # plot without X
        plt.scatter(np.ravel(hqcs_no_x_freqs), np.ravel(ccs_no_x_freqs), alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('CCS frequency')
        plt.title("{} no X".format(label))
        plt.savefig(out("freq_agreement_no_x_{}.pdf".format(label)))
        plt.close()

        js_divs = np.array(list(js_divergence(a, b)
                                for a, b in zip(hqcs_no_x_freqs.T, ccs_no_x_freqs.T)))

        bad_columns = np.where(js_divs > cutoff)[0]

        # plot column-wise JS divergence
        plt.figure(figsize=(8, 5))
        plt.plot(js_divs)
        if len(bad_columns) > 0:
            plt.axhline(y=cutoff, color='r')
        plt.xlabel('column')
        plt.ylabel('JS divergence')
        plt.savefig(out("js_div_vals_{}.pdf".format(label)), )
        plt.close()

        # plot histogram of JS divergence
        plt.figure(figsize=(8, 5))
        plt.hist(js_divs)
        plt.xlabel('JS divergence')
        plt.ylabel('number of columns')
        plt.savefig(out("js_div_hist_{}.pdf".format(label)))
        plt.close()

        # compute JS divergence and frequence of top base
        cols = [np.arange(len(js_divs)),
                js_divs,
                keys[hqcs_no_x_freqs.argmax(axis=0)],
                hqcs_no_x_freqs.max(axis=0),
                keys[ccs_no_x_freqs.argmax(axis=0)],
                ccs_no_x_freqs.max(axis=0)]
        rows = list(zip(*cols))
        with open(out("js_divergence_{}.csv".format(label)), "w") as f:
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
