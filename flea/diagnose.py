#!/usr/bin/env python

"""
Compare HQCS frequencies to QCS frequencies.

Usage:
  diagnose.py [options] <hqcs_file> <qcs_file> <result_path>
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

from flea.util import column_count, prob, js_divergence
from flea.util import id_to_label, id_to_cn


np.set_printoptions(suppress=True)


def diagnose(hqcsfile, qcsfile, result_path, cutoff):
    def out(f):
        return os.path.join(result_path, f)

    format = "fasta"

    plt.style.use('seaborn-deep')

    #Getting and sorting HQCS seqs
    hqcs = list(sorted(SeqIO.parse(hqcsfile, format), key=lambda r: r.id))
    hqcs_array = np.array(list(list(str(rec.seq)) for rec in hqcs))
    hqcs_labels = np.array(list(id_to_label(rec.id) for rec in hqcs))
    copynumber_array = np.array(list(id_to_cn(rec.id) for rec in hqcs))

    # import qcs seqs
    qcs = list(SeqIO.parse(qcsfile, format))
    qcs_array = np.array(list(list(str(rec.seq)) for rec in qcs))
    qcs_labels = np.array(list(id_to_label(rec.id) for rec in qcs))

    if hqcs_array.shape[1] != qcs_array.shape[1]:
        raise Exception('alignment shapes do not match')

    # amino acids
    aas = list(sorted(set(hqcs_array.ravel()) | set(qcs_array.ravel())))
    keys = np.array(aas)

    # Relying on "_" to separate timepoint ID
    labels = sorted(set(hqcs_labels))

    for label in labels + ['overall']:
        if label == 'overall':
            hqcs_bools = np.array([True] * len(hqcs_array))
            qcs_bools = np.array([True] * len(qcs_array))
        else:
            hqcs_bools = hqcs_labels == label
            qcs_bools = qcs_labels == label

        hqcs_counts = column_count(hqcs_array[hqcs_bools], keys, weights=copynumber_array[hqcs_bools])
        qcs_counts = column_count(qcs_array[qcs_bools], keys)

        if np.all(hqcs_counts.sum(axis=0) != copynumber_array[hqcs_bools].sum()):
            raise Exception('bad hqcs counts')

        if np.all(qcs_counts.sum(axis=0) != qcs_bools.sum()):
            raise Exception('bad qcs counts')

        np.savetxt(out("hqcs_counts_{}.csv".format(label)), hqcs_counts.T,
                   fmt="%1.1u", delimiter=",", header=",".join(aas))
        np.savetxt(out("qcs_counts_{}.csv".format(label)), qcs_counts.T,
                   fmt="%u", delimiter=",", header=",".join(aas))

        # plot all freqs
        plt.scatter(np.ravel(prob(hqcs_counts, axis=0)),
                    np.ravel(prob(qcs_counts, axis=0)),
                    alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('QCS frequency')
        plt.title(label)
        plt.savefig(out("freq_agreement_{}.pdf".format(label)))
        plt.close()

        # strip X and normalize
        x_index = aas.index('X')
        hqcs_no_x = np.delete(hqcs_counts, x_index, 0)
        hqcs_no_x_freqs = prob(hqcs_no_x, axis=0)
        qcs_no_x = np.delete(qcs_counts, x_index, 0)
        qcs_no_x_freqs = prob(qcs_no_x, axis=0)

        # plot without X
        plt.scatter(np.ravel(hqcs_no_x_freqs), np.ravel(qcs_no_x_freqs), alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('QCS frequency')
        plt.title("{} no X".format(label))
        plt.savefig(out("freq_agreement_no_x_{}.pdf".format(label)))
        plt.close()

        js_divs = np.array(list(js_divergence(a, b)
                                for a, b in zip(hqcs_no_x_freqs.T, qcs_no_x_freqs.T)))

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
                keys[qcs_no_x_freqs.argmax(axis=0)],
                qcs_no_x_freqs.max(axis=0)]
        rows = list(zip(*cols))
        with open(out("js_divergence_{}.csv".format(label)), "w") as f:
            writer = csv.writer(f)
            writer.writerow(['column', 'js_divergence', 'hqcs_aa', 'hqcs_freq', 'qcs_aa', 'qcs_freq'])
            writer.writerows(rows)


if __name__ == "__main__":
    args = docopt(__doc__)
    hqcs_file = args["<hqcs_file>"]
    qcs_file = args["<qcs_file>"]
    result_path = args["<result_path>"]
    cutoff = float(args["-c"])
    diagnose(hqcs_file, qcs_file, result_path, cutoff)
