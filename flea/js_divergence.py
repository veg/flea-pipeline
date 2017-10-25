import sys
import csv
import json

import click
from Bio import SeqIO
import numpy as np

from flea.util import js_divergence
from flea.util import column_count
from flea.util import prob
from flea.util import get_date_dict
from flea.util import id_to_copynumber, id_to_label
from flea.util import read_single_record


@click.command()
@click.argument('infile')
@click.argument('mrcafile')
@click.argument('metafile')
@click.argument('outfile')
def main(infile, mrcafile, metafile, outfile):
    date_dict = get_date_dict(metafile)

    records = list(SeqIO.parse(infile, 'fasta'))

    seqs = list(list(str(r.seq)) for r in records)
    # assume id has "<timepoint>_<misc>_n_<copynumber>" pattern
    copynumbers = list(id_to_copynumber(r.id) for r in records)
    seq_dates = list(date_dict[id_to_label(r.id)] for r in records)

    # add MRCA
    mrca_record = read_single_record(mrcafile, 'fasta', expected=1)
    seqs.insert(0, list(str(mrca_record.seq)))
    copynumbers.insert(0, 1)
    seq_dates.insert(0, 0)

    seq_array = np.array(seqs)
    copynumber_array = np.array(copynumbers)
    date_array = np.array(seq_dates)

    # FIXME: labels must sort according to date!!!
    dates = list(sorted(set(date_dict.values())))
    dates.insert(0, 0)

    aas = list(sorted(set(seq_array.ravel())))
    alphabet = np.array(aas)

    ps = {}
    for date in dates:
        bools = (date_array == date)
        counts = column_count(seq_array[bools], alphabet,
                              weights=copynumber_array[bools])
        ps[date] = prob(counts, axis=0).T

    result = {}
    for i in range(1, len(dates)):
        prev_date = dates[i - 1]
        cur_date = dates[i]
        prev_ps = ps[prev_date]
        cur_ps = ps[cur_date]
        result[dates[i]] = list(js_divergence(a, b)
                                for a, b in zip(prev_ps, cur_ps))
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",", ":"))


if __name__ == "__main__":
    main()
