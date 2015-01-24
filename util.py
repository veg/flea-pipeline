"""Utility functions used in more than one script."""

from subprocess import Popen, PIPE, STDOUT
from itertools import zip_longest
from itertools import tee
from itertools import filterfalse
import time
import os

from Bio.Seq import Seq


def call(cmd_str, stdin=None, stdout=None):
    """Call a command in the shell. Captures STDOUT and STDERR.

    If *args are given, they are passed as strings to STDIN.

    Raises an exception if return code != 0.

    """
    if isinstance(stdin, str):
        in_str = stdin.encode()
        stdin = PIPE
    else:
        in_str = None
    if stdout is None:
        stdout = PIPE
    process = Popen(cmd_str, stdin=stdin, stdout=stdout, stderr=PIPE, shell=True)
    stdout_str, stderr_str = process.communicate(input=in_str)
    if process.returncode != 0:
        if in_str is None:
            in_str = ''
        else:
            in_str = in_str.decode()
        raise Exception("Failed to run '{}'\nExit status: {}\nSTDIN:\n{}\nSTDOUT:\n{}\nSTDERR\n{}\n".format(
                cmd_str, process.returncode, in_str, stdout_str.decode(), stderr_str.decode()))


def qsub(cmd, sentinel, walltime=3600, sleep=1):
    """A blocking qsub.

    sentinel: file to be created when task is done. Cannot exist.
    walltime: seconds

    """
    if walltime < 1:
        raise Exception('walltime={} < 1'.format(walltime))
    if os.path.exists(sentinel):
        raise Exception('sentinel already exists!')
    mycmd = '{}; echo $? > {}'.format(cmd, sentinel)
    fwalltime = time.strftime('%H:%M:%S', time.gmtime(walltime))
    qsub_cmd = 'qsub -V -W umask=077 -l walltime={}'.format(fwalltime)
    full_cmd = 'echo "{}" | {}'.format(mycmd, qsub_cmd)
    call(full_cmd)
    run_time = 0
    while run_time < walltime:
        time.sleep(sleep)
        run_time += sleep
        if os.path.exists(sentinel):
            break
    # wait a second to make sure it has been flushed
    # TODO: this is probably not robust
    time.sleep(1)
    with open(sentinel) as handle:
        code = handle.read().strip()
        if code != '0':
            raise Exception('qsub job "{}" exited with code "{}"'.format(full_cmd, code))


def hyphy_call(script_file, *args, hyphy=None):
    if hyphy is None:
        hyphy = "HYPHYMP"
    cmd = '{} {}'.format(hyphy, script_file)
    if args:
        in_str = "".join("{}\n".format(i) for i in args)
    else:
        in_str = None
    call(cmd, stdin=in_str)


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


def flatten(it):
    return (b for a in it for b in a)


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


def touch(f):
    call('touch {}'.format(f))


def strlist(arg):
    if isinstance(arg, str):
        return [arg]
    return arg


def write_to_handle(handle, output):
    if isinstance(handle, str):
        handle = open(handle, 'w')
        do_close = True
    try:
        handle.write(output)
    finally:
        if do_close:
            handle.close()
