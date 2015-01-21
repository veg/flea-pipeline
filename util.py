"""Utility functions used in more than one script."""

from itertools import zip_longest

from Bio.Seq import Seq


def insert_gaps(source, target, src_gap, target_gap, skip):
    """Inserts `target_gap` into target string at same positions as
    `src_gap` in `source`.

    >>> insert_gaps("a-bc", "def", "-", "-", 1)
    'd-ef'

    >>> insert_gaps("a-bc", "def", "-", "--", 1)
    'd--ef'

    >>> insert_gaps("a-bc", "dddeeefff", "-", "---", 3)
    'ddd---eeefff'

    """
    if len(src_gap) != 1:
        # TODO: remove this restriction
        raise ValueError("Argument `src_gap` must be a single character."
                         " Current value: '{}'".format(src_gap))
    if len(list(c for c in source if c != src_gap)) * skip != len(target):
        raise ValueError("`src` and `target` have different lengths.")
    result = []
    for c in source:
        if c == src_gap:
            result.append(target_gap)
        else:
            result.append(target[0:skip])
            target = target[skip:]
    return "".join(result)
    

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def new_record_seq_str(record, seqstr):
    record = record[:]
    record.seq = Seq(seqstr, alphabet=record.seq.alphabet)
    return record
