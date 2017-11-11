"""Utility functions used in more than one script."""

from fnmatch import fnmatch
from itertools import zip_longest
from itertools import tee
from itertools import filterfalse
from itertools import islice
import re
import time
import datetime
import os
from functools import wraps
from glob import glob
import sys
from contextlib import contextmanager
import csv
import random

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Seq import translate

import numpy as np


def insert_gaps(source, target, src_gap, target_gap, gap_char=None):
    """Inserts `target_gap` into target string at same positions as
    `src_gap` in `source`.

    >>> insert_gaps("a-bc", "def", "-", "-")
    'd-ef'

    >>> insert_gaps("a-bc", "ddeeff", "-", "--")
    'dd--eeff'

    >>> insert_gaps("a-bc", "dddeeefff", "-", "---")
    'ddd---eeefff'

    >>> insert_gaps("a--bc", "dddeeefff", "-", "---")
    'ddd------eeefff'

    >>> insert_gaps("aaa---bbbccc", "abc", "---", "-")
    'a-bc'

    """
    if gap_char is None:
        gap_char = '-'
    if not all(g == gap_char for g in src_gap):
        raise ValueError('illegal source gap')
    if not all(g == gap_char for g in target_gap):
        raise ValueError('illegal target gap')
    src_skip = len(src_gap)
    target_skip = len(target_gap)
    if len(source) % src_skip != 0:
        raise ValueError('source length does not match gap length')
    if len(target) % target_skip != 0:
        raise ValueError('target length does not match gap length')

    ungapped_src = ''.join(list(''.join(chunk) for chunk in grouper(source, src_skip)
                                if ''.join(chunk) != src_gap))
    if len(ungapped_src) % src_skip != 0:
        raise ValueError('ungapped source length does not match gap length')
    if len(ungapped_src) // src_skip != len(target) // target_skip:
        raise ValueError("`src` and `target` have different lengths.")
    result = []
    for chunk in grouper(source, src_skip):
        c = ''.join(chunk)
        if c == src_gap:
            result.append(target_gap)
        elif gap_char in chunk:
            raise ValueError('{} contains a gap character'.format(chunk))
        else:
            result.append(target[0:target_skip])
            target = target[target_skip:]
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


def nth(iterable, n, default=None):
    "Returns the nth item or a default value"
    return next(islice(iterable, n, None), default)


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
    record.letter_annotations = {}
    record.seq = Seq(seqstr, alphabet=record.seq.alphabet)
    return record


def update_record_seq(record, seq):
    record.seq = seq
    return record


def flatten(it):
    return (b for a in it for b in a)


def traverse(o, tree_types=(list, tuple)):
    if isinstance(o, tree_types):
        for value in o:
            for subvalue in traverse(value):
                yield subvalue
    else:
        yield o


def cat_files(files, outfile, chunk=2**14):
    """Concatenate multiple files in chunks."""
    out_handle = open(outfile, "w")
    for f in files:
        handle = open(f)
        while True:
            data = handle.read(chunk)
            if data:
                out_handle.write(data)
            else:
                break


def strlist(arg):
    if isinstance(arg, str):
        return [arg]
    return arg


def ensure_not_empty(files):
    for f in strlist(files):
        if not os.path.exists(f):
            raise Exception('Expected output file does not'
                            ' exist: "{}"'.format(f))
        if os.stat(f).st_size == 0:
            raise Exception('Empty file: "{}"'.format(f))


def read_single_record(filename, format, expected=None):
    records = list(SeqIO.parse(filename, format))
    if expected is not None and len(records) != expected:
        raise Exception('expected only {} records in {}'.format(expected, filename))
    return records[0]


def parse_by_ext(filename):
    format = filename[-len('fastx'):]
    return SeqIO.parse(filename, format)


def check_seq_number(filename, min_n):
    found_n = sum(1 for _ in parse_by_ext(filename))
    if found_n < min_n:
        raise Exception('file "{}" must have at least {} entries,'
                        'but it only has {}'.format(filename, min_n, found_n))


def check_seq_ids(inputs, output):
    a_ids = set(r.id for a in inputs for r in parse_by_ext(a))
    b_ids = set(r.id for r in parse_by_ext(output))
    unmatched = a_ids.symmetric_difference(b_ids)
    if unmatched:
        raise Exception('IDs from "{}" do not match "{}".'
                        ' Number of unmatched ids: {}'.format(inputs, output, len(unmatched)))


def check_unique_ids(f):
    ids = list(r.id for r in parse_by_ext(f))
    if len(set(ids)) < len(ids):
        raise Exception('IDs in "{}" are not'
                        ' unique'.format(f))


def check_seq_ratio(inputs, output, expected):
    a_exp, b_exp = expected
    if a_exp != int(a_exp):
        raise Exception('expected an int, but got {}'.format(a_exp))
    if b_exp != int(b_exp):
        raise Exception('expected an int, but got {}'.format(b_exp))
    a_n = sum(1 for a in inputs for r in parse_by_ext(a))
    b_n = sum(1 for r in parse_by_ext(output))
    if a_n * b_exp != b_n * a_exp:
        raise Exception('Sequnce ratios do not match.'
                        ' Expected {}:{}. Got: {}:{}.'
                        ' Inputs: "{}" Output: "{}"'.format(
                a_exp, b_exp, a_n, b_n, inputs, output))


def check_illegal_chars(f, chars):
    for r in parse_by_ext(f):
        found = set(chars.upper()).intersection(set(str(r.seq).upper()))
        if found:
            raise Exception('Illegal characters "{}" found in sequence "{}"'
                            ' of file "{}"'.format(found, r.id, f))


def check_in_frame(f):
    for r in parse_by_ext(f):
        len_ = len(r.seq.ungap('-'))
        if len_ % 3 != 0:
            raise Exception('Sequence "{}" in file "{}" has length={} not multiple'
                            ' of 3'.format(r.id, f, len_))


def extend_coordinates(coordinates, seq, gap=None):
    """Extend coordinates to a gappy sequence.

    >>> extend_coordinates([1, 2, 3, 4], "a-b-cd")
    [1, 1, 2, 2, 3, 4]

    """
    if gap is None:
        gap = "-"
    if sum(1 for c in seq if c != gap) != len(coordinates):
        raise Exception('coordinates do not match source')
    coords_gen = iter(coordinates)
    coord = None
    result = []
    for char in seq:
        if char != gap:
            coord = next(coords_gen)
        result.append(1 if coord is None else coord)
    return result


def column_count(a, alphabet, weights=None):
    """
    >>> column_count(np.array([['A', 'A'], ['A', 'B']]), ['A', 'B'])
    array([[2, 1],
           [0, 1]])

    >>> column_count(np.array([['A', 'A'], ['A', 'B']]), alphabet=['A', 'B'], weights=[3, 2])
    array([[5, 3],
           [0, 2]])

    """
    alphabet = np.array(alphabet).ravel()
    result = (a == alphabet[:, np.newaxis, np.newaxis])
    if weights is not None:
        weights = np.array(weights).ravel()
        assert len(weights) == len(a)
        result = result * weights.reshape(-1, 1)
        assert np.all(result.sum(axis=1).sum(axis=0) == weights.sum())
    return result.sum(axis=1)


def prob(p, axis=None):
    result = np.array(p)
    if result.sum() != 1:
        result = result / result.sum(axis=axis)
    return result


def entropy(p):
    """
    >>> np.abs(entropy([0, 1]))
    0.0

    >>> entropy([0.5, 0.5])
    1.0

    >>> entropy([0.25, 0.25, 0.25, 0.25])
    2.0

    """
    p = prob(p)
    p = p[p > 0]
    result = -(p * np.log2(p)).sum()
    if np.any(np.isnan(result)):
        raise Exception('entropy failed')
    return result


def js_divergence(*args, weights=None):
    """Jensen-Shannon divergence

    >>> js_divergence([0.25, 0.25, 0.5], [0.25, 0.25, 0.5])
    0.0

    >>> js_divergence([1.0, 0.0], [0.0, 1.0])
    1.0


    """
    # ensure there are no zeros
    if weights is None:
        weights = prob(np.repeat(1, len(args)))
    ps = np.vstack(list(prob(p) for p in args))
    m = (weights.reshape(-1, 1) * ps).sum(axis=0)
    bools = (m > 0)
    ps = ps[:, bools]
    m = m[bools]
    result = entropy(m) - (weights * np.apply_along_axis(entropy, axis=1, arr=ps)).sum()
    if result < 0:
        raise Exception('JS divergence cannot be 0')
    return np.abs(result)  # to avoid -0


def str_to_type(x):
    """Try to convert string to Python types.

    There are probably more general ways to do this, but ints, floats,
    and bools are good enough for our purposes.

    >>> str_to_type("3.0")
    3.0

    >> str_to_type("3")
    3

    >> str_to_type("True")
    True

    >> str_to_type("false")
    False

    >> str_to_type("a_string")
    "a_string"

    """
    try:
        if str(int(x)) == x:
            return int(x)
    except ValueError:
        pass
    try:
        return float(x)
    except ValueError:
        pass
    if x.lower() == "true":
        return True
    if x.lower() == "false":
        return False
    return x


def run_regexp(runlen, targets=None):
    """Get compiled regexp for `runlen` runs of characters.

    >>> bool(run_regexp(3).search("AAA"))
    True

    >>> bool(run_regexp(3).search("AATT"))
    False

    >>> bool(run_regexp(3, "T").search("AAATT"))
    False

    """
    if targets is None:
        targets = "ACGT"
    pattern = "|".join("({target}){{{runlen}}}".format(target=t, runlen=runlen)
                       for t in targets)
    return re.compile(pattern)


def iter_sample(iterable, k):
    """Sample `k` items from an iterable.

    Adapted from http://stackoverflow.com/a/12581484

    """
    results = []
    for i, v in enumerate(iterable):
        r = random.randint(0, i)
        if r < k:
            if i < k:
                results.insert(r, v) # add first `n` items in random order
            else:
                results[r] = v # at a decreasing rate, replace random items
    if len(results) < k:
        raise ValueError("Sample larger than population.")
    return results


def get_format(f):
    return os.path.splitext(f)[1][1:]


def is_in_frame(seq, allow_stop_codons):
    if len(seq) % 3 != 0:
        return False
    t = seq.translate()
    if len(t) != len(seq) / 3:
        return False
    if '*' in t and not allow_stop_codons:
        return False
    return True


def parse_copynumbers(infile):
    with open(infile) as handle:
        parsed = csv.reader(handle, delimiter='\t')
        return dict((i, int(n)) for i, n in parsed)


def replace_id(record, id_):
    result = record[:]
    result.id = id_
    # HyPhy fails if a name or description are present
    result.name = ""
    result.description = ""
    return result


def id_with_copynumber(id_, copynumber):
    return "{}_n_{}".format(id_, copynumber)


def id_to_copynumber(id_):
    return int(id_.split('_')[-1])


def id_to_label(id_):
    return id_.split('_')[0]


def get_date_dict(metafile):
    with open(metafile) as handle:
        parsed = csv.reader(handle, delimiter=' ')
        return dict((tp, int(date)) for _, tp, date in parsed)


def handle_codon(codon):
    codon = ''.join(codon)
    if codon == "---":
        return codon
    if translate(codon) == "*":
        return "NNN"
    return codon


def replace_stop_codons(record):
    if len(record.seq) % 3 != 0:
        raise Exception('record {} is not in frame'.format(record.id))
    new_str = "".join(handle_codon(codon) for codon in grouper(record.seq, 3))
    result = new_record_seq_str(record, new_str)
    # HyPhy does not like anything but the sequence id
    result.name = ""
    result.description = ""
    return result


def handle_gap_codon(codon):
    codon = ''.join(codon)
    if '-' in codon and codon != '---':
        return 'NNN'
    return codon


def replace_gapped_codons(record):
    if len(record.seq) % 3 != 0:
        raise Exception('record {} is not in frame'.format(record.id))
    new_str = "".join(handle_gap_codon(codon) for codon in grouper(record.seq, 3))
    result = new_record_seq_str(record, new_str)
    # HyPhy does not like anything but the sequence id
    result.name = ""
    result.description = ""
    return result

