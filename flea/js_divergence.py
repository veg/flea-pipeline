import sys
import csv
import json

from Bio import SeqIO
import numpy as np

from flea.util import js_divergence
from flea.util import column_count
from flea.util import prob
from flea.util import get_date_dict


def main(infile, metafile, outfile):
    date_dict = get_date_dict(metafile)
    
    records = list(SeqIO.parse(infile, 'fasta'))
    seq_array = np.array(list(list(str(r.seq)) for r in records))
    # assume id has "<timepoint>_<misc>_cn_<copynumber>" pattern
    copynumber_array = np.array(list(int(r.id.split("_")[-1])
                                     for r in records))
    date_array = np.array(list(date_dict[r.id.split("_")[0]]
                               for r in records))

    # FIXME: labels must sort according to date!!!
    dates = list(sorted(set(date_dict.values())))

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
    infile = sys.argv[1]
    metafile = sys.argv[2]
    outfile = sys.argv[3]
    main(infile, metafile, outfile)
