"""Utility functions used in more than one script."""

from itertools import zip_longest
from itertools import tee
from itertools import filterfalse

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


def partition(pred, iterable):
    'Use a predicate to partition entries into false entries and true entries'
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = tee(iterable)
    return filterfalse(pred, t1), filter(pred, t2)


def genlen(gen, ldict, name):
    """A bit of a hack to get the length of a generator after its been run.

    ldict = {}
    list(genlen((i for i in range(3)), ldict, 'mylen')
    assert(ldict['mylen'] == 3

    """
    ldict[name] = 0
    def f():
        while True:
            yield next(gen)
            ldict[name] += 1
    return f()


def new_record_seq_str(record, seqstr):
    record = record[:]
    record.seq = Seq(seqstr, alphabet=record.seq.alphabet)
    return record
