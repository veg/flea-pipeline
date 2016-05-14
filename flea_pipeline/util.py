"""Utility functions used in more than one script."""

from fnmatch import fnmatch
from subprocess import Popen, PIPE, CalledProcessError
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
import warnings
import sys

from Bio.Seq import Seq
from Bio import SeqIO

import numpy as np

from ruffus.drmaa_wrapper import run_job, error_drmaa_job

import flea_pipeline.pipeline_globals as globals_


def format_walltime(seconds):
    seconds = int(seconds)
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{:02}:{:02}:{:02}".format(h, m, s)


# TODO: rename
# TODO: handle stdout and stderr
def maybe_qsub(cmd, infiles, outfiles, ppn=1, walltime=None,
               stdout=None, stderr=None, name=None):
    from uuid import uuid4  # because this hangs on silverback compute nodes
    if walltime is None:
        walltime = globals_.config.getint('Jobs', 'walltime')
    if name is None:
        name = 'job-{}'.format(os.path.basename(cmd.split()[0]))
    name = "{}-{}".format(name, uuid4())

    queue_name = globals_.config.get('Jobs', 'queue')
    fwalltime = format_walltime(walltime)
    nodes = 1
    job_other_options = ('-q {} -l walltime={},nodes={}:ppn={}'.format(queue_name, fwalltime, nodes, ppn))

    success = False
    msg = ''
    stdout_res, stderr_res = "", ""
    try:
        stdout_res, stderr_res = run_job(cmd_str=cmd,
                                         job_name=name,
                                         logger=globals_.logger,
                                         drmaa_session = globals_.drmaa_session,
                                         run_locally=globals_.run_locally,
                                         job_other_options=job_other_options,
                                         job_script_directory=globals_.job_script_dir)
        success = True
    except error_drmaa_job as err:
        msg = err.msg

    return Result(infiles, outfiles, stdout=stdout_res, stderr=stderr_res,
                  desc=name, failed=(not success),
                  msg=msg, cmds=cmd)



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
                        ' Unmatched ids: {}'.format(inputs, output, sorted(unmatched)))


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

class Result(object):
    # borrowed from https://github.com/daler/pipeline-example
    def __init__(self, infiles, outfiles, log=None, stdout=None, stderr=None,
                 desc=None, failed=False, msg=None, cmds=None):
        """
        A Result object encapsulates the information going to and from a task.
        Each task is responsible for determining the following arguments, based
        on the idiosyncracies of that particular task:
            `infiles`: Required list of input files to the task
            `outfiles`: Required list of output files to the task
            `failed`: Optional (but recommended) boolean indicating whether the
                      task failed or not.  Default is False, implying that the
                      task will always work
            `log`: Optional log file from running the task's commands
            `stdout`: Optional stdout that will be printed in the report if the
                      task failed
            `stderr`: Optional stderr that will be printed in the report if the
                      task failed
            `desc`: Optional description of the task
            `cmds`: Optional string of commands used by the task. Can be very
                    useful for debugging.
        """
        if isinstance(infiles, str):
            infiles = [infiles]
        if isinstance(outfiles, str):
            outfiles = [outfiles]
        self.infiles = infiles
        self.outfiles = outfiles
        self.log = log
        self.stdout = stdout
        self.stderr = stderr
        self.elapsed = None
        self.failed = failed
        self.msg = msg
        self.desc = desc
        self.cmds = cmds

    def report(self, logger_proxy, logging_mutex):
        """Prints a nice report."""
        with logging_mutex:
            if not self.desc:
                self.desc = ""
            logger_proxy.info('Finish task: {}'.format(self.desc))
            logger_proxy.info('    Time: {}'.format(datetime.datetime.now()))
            if self.elapsed is not None:
                logger_proxy.info('    Elapsed: {}'.format(nicetime(self.elapsed)))
            if self.cmds is not None:
                logger_proxy.info('    Commands: {}'.format(str(self.cmds)))
            for input_fn in self.infiles:
                input_fn = os.path.normpath(os.path.relpath(input_fn))
                logger_proxy.info('    Input: {}'.format(input_fn))
            for output_fn in self.outfiles:
                output_fn = os.path.normpath(os.path.relpath(output_fn))
                logger_proxy.info('    Output: {}'.format(output_fn))
            if self.log is not None:
                logger_proxy.info('    Log: {}'.format(self.log))
            if self.failed:
                logger_proxy.error('=' * 80)
                logger_proxy.error('Error in {}'.format(self.desc))
                if self.msg:
                    logger_proxy.error('message: {}'.format(self.msg))
                if self.log is not None:
                    logger_proxy.error('log: {}'.format(self.log))
                if self.stderr:
                    logger_proxy.error('====STDERR====')
                    logger_proxy.error(self.stderr)
                if self.stdout:
                    logger_proxy.error('====STDOUT====')
                    logger_proxy.error(self.stdout)
                logger_proxy.error('=' * 80)
            logger_proxy.info('')
            if self.failed:
                sys.exit(1)



def nicetime(seconds):
    """Convert seconds to hours-minutes-seconds"""
    # borrowed from https://github.com/daler/pipeline-example
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    elapsed = "%dh%02dm%02.2fs" % (h, m, s)
    return elapsed


def report_wrapper(function):
    @wraps(function)  # necessary because ruffus uses function name internally
    def wrapped(infiles, outfiles, *args, **kwargs):
        # TODO: log start, infiles, outfiles, cmd
        t0 = time.time()
        result = function(infiles, outfiles, *args, **kwargs)
        t1 = time.time()
        if result is None:
            result = Result(infiles=infiles, outfiles=outfiles)
        result.elapsed = (t1 - t0)
        if not result.desc:
            result.desc = function.__name__
        result.report(globals_.logger, globals_.logger_mutex)
    return wrapped


def must_work(maybe=False, seq_ratio=None, seq_ids=False, min_seqs=None,
              unique_ids=False, many=False,
              illegal_chars=None, pattern=None, in_frame=False):
    """Fail if any output is empty.

    seq_ratio: (int, int): ensure numbers of sequences match
    seq_ids: ensure all ids present in input files are in output file
    min_seqs: minimum number of sequences in the output file
    unique_ids: sequence ids must be unique
    many: produces an unknown number of output files. Search with a pattern.
    pattern: pattern for determining input files to consider

    """
    if pattern is None:
        pattern = '*'

    def wrap(function):
        if globals_.options.touch_files_only or globals_.options.just_print:
            return function

        @wraps(function)  # necessary because ruffus uses function name internally
        def wrapped(infiles, outfiles, *args, **kwargs):
            function(infiles, outfiles, *args, **kwargs)
            if seq_ids or seq_ratio:
                infiles = list(traverse(strlist(infiles)))
                infiles = list(f for f in infiles if fnmatch(f, pattern))
            if many:
                outdir, outpattern = args[-2:]
                outfiles = list(os.path.join(outdir, f)
                                for f in os.listdir(outdir) if re.search(outpattern, f))
                if not outfiles:
                    raise Exception('Task produced no output')
            outfiles = strlist(outfiles)
            ensure_not_empty(outfiles)
            if min_seqs is not None:
                check_seq_number(outfiles[0], min_seqs)
            if seq_ids:
                if len(outfiles) != 1:
                    raise Exception('cannot check seq_ids for'
                                    ' more than one outfile')
                check_seq_ids(infiles, outfiles[0])
            if unique_ids:
                for o in outfiles:
                    check_unique_ids(o)
            if seq_ratio is not None:
                if len(outfiles) != 1:
                    raise Exception('cannot check seq_ratio for'
                                    ' more than one outfile')
                check_seq_ratio(infiles, outfiles[0], seq_ratio)
            if illegal_chars:
                for f in outfiles:
                    check_illegal_chars(f, illegal_chars)
            if in_frame:
                for f in outfiles:
                    check_in_frame(f)
        return wrapped
    return wrap


def must_produce(n):
    """Check that at least `n` files match the pathname glob"""
    def wrap(function):
        @wraps(function)
        def wrapped(infiles, outfiles, pathname, *args, **kwargs):
            function(infiles, outfiles, pathname, *args, **kwargs)
            n_produced = len(glob(pathname))
            if n_produced < n:
                raise Exception('Task was supposed to produce at least {}'
                                ' outputs, but it produced {}'.format(n, n_produced))
        return wrapped
    return wrap


def check_suffix(name, suffix):
    if not name.endswith(suffix):
        raise Exception("filename '{}' does not "
                        " end with '{}'".format(name, suffix))


def remove_suffix(name, suffix):
    check_suffix(name, suffix)
    return name[:-len(suffix)]


def check_basename(name, bn):
    exp_bn = os.path.basename(name)
    if exp_bn != bn:
        raise Exception("filename '{}' does not"
                        " match expected name '{}'".format(bn, exp_bn))


def split_name(name):
    # TODO: duplicate code with name_to_label
    cands = list(t.key for t in globals_.timepoints
                 if name.startswith(t.key))
    if not cands:
        raise Exception('name matches no'
                        ' key: "{}"'.format(name))
    if len(cands) != 1:
        raise Exception('name starts with non-unique'
                        ' key: "{}"'.format(name))
    key = cands[0]
    rest = name[len(key):].strip('_')
    if not rest:
        raise Exception('name has no unique part:'
                        ' "{}"'.format(name))
    return key, rest


def name_key_to_label(name):
    key, _ = split_name(name)
    return globals_.key_to_label[key]


def name_to_label(name):
    cands = list(t.label for t in globals_.timepoints
                 if name.startswith(t.label))
    if not cands:
        raise Exception('name matches no'
                        ' label: "{}"'.format(name))
    if len(cands) != 1:
        raise Exception('name starts with non-unique'
                        ' label: "{}"'.format(name))
    label = cands[0]
    rest = name[len(label):].strip('_')
    if not rest:
        raise Exception('name has no unique part:'
                        ' "{}"'.format(name))
    return label


def name_to_date(name):
    return globals_.label_to_date[name_to_label(name)]


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
    cur = -1
    result = []
    for char in seq:
        if char != gap:
            cur = next(coords_gen)
        result.append(cur)
    return result


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


def translate_helper(infile, outfile, gapped, name):
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "translate.py"),
        'infile': infile,
        'outfile': outfile,
        'gap_option': '--gapped' if gapped else '',
        }
    cmd = "{python} {script} {gap_option} {infile} {outfile}".format(**kwargs)
    return maybe_qsub(cmd, infile, outfile, name=name)


def usearch_hqcs_ids(infile, outfile, dbfile, name=None):
    """run usearch_global against a hqcs database and print pairs of ids"""
    identity = globals_.config.get('Parameters', 'ccs_to_hqcs_identity')
    max_accepts = globals_.config.get('Parameters', 'max_accepts')
    max_rejects = globals_.config.get('Parameters', 'max_rejects')
    maxqt = globals_.config.get('Parameters', 'max_query_target_length_ratio')
    cmd = ("{usearch} -usearch_global {infile} -db {db} -id {id}"
           " -userout {outfile} -top_hit_only -userfields query+target -strand both"
           " -maxaccepts {max_accepts} -maxrejects {max_rejects}"
           " -maxqt {maxqt}")
    if name is None:
        name = 'usearch-hqcs-ids'
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, db=dbfile, id=identity, outfile=outfile,
                     max_accepts=max_accepts, max_rejects=max_rejects,
                     maxqt=maxqt)
    return maybe_qsub(cmd, infile, outfiles=outfile, name=name)
