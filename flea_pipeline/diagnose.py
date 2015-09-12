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
        plt.scatter(np.ravel(prob(hqcs_counts, axis=0)),
                    np.ravel(prob(ccs_counts, axis=0)),
                    alpha=0.5)
        plt.xlabel('HQCS frequency')
        plt.ylabel('CCS frequency')
        plt.title(label)
        plt.savefig(out("freq_agreement_{}.png".format(label)))
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
