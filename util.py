"""Utility functions used in more than one script."""

from subprocess import Popen, PIPE, STDOUT
from itertools import zip_longest
from itertools import tee
from itertools import filterfalse
from collections import defaultdict
from uuid import uuid4
import time
import os

from Bio.Seq import Seq


# TODO: capture and write all stdout and stderr to files.
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
        raise Exception("Failed to run '{}'\nExit status: {}\nSTDIN:\n{}"
                        "\nSTDOUT:\n{}\nSTDERR\n{}\n".format(
                cmd_str, process.returncode, in_str, stdout_str.decode(), stderr_str.decode()))
    return stdout_str.decode()


def job_state(job_id):
    """Returns job state from qstat.

    If the job is not found, returns 'M' for 'missing'.

    """
    cmd = 'qstat {}'.format(job_id)
    try:
        result = call(cmd)
        name, _, _, _, state, _ = result.split('\n')[2].split()
        return state
    except CalledProcessError:
        pass
    return 'M'


def wait_for_job(job_id, sleep):
    """Check job state every `sleep` seconds; return if 'C' or 'M'."""
    while True:
        time.sleep(sleep)
        state = job_state(job_id)
        if state in 'CM':
            return


def wait_for_files(files, sleep, walltime):
    """Wait up to `walltime` seconds for files in `files` to exist"""
    files = strlist(files)
    run_time = 0
    while run_time < walltime:
        time.sleep(sleep)
        run_time += sleep
        if all(os.path.exists(f) for f in files):
            break


def format_walltime(seconds):
    seconds = int(seconds)
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{:02}:{:02}:{:02}".format(h, m, s)


def qsub(cmd, outfiles=None, queue=None, nodes=1, ppn=1, sleep=5,
         walltime=3600, waittime=10, name=None, stdout=None, stderr=None):
    """A blocking qsub."""
    if walltime < 1:
        raise Exception('walltime={} < 1'.format(walltime))
    if name is None:
        name = 'job-{}'.format(os.path.basename(cmd.split()[0]))
    sentinel = '{}-{}.complete'.format(name, uuid4())
    mycmd = '{}; echo \$? > {}'.format(cmd, sentinel)
    fwalltime = format_walltime(walltime)
    qsub_cmd = ('qsub -V -W umask=077 -l'
                ' walltime={},nodes={}:ppn={} -N {name}'
                ' -d `pwd` -w `pwd`'.format(fwalltime, nodes, ppn, name=name))
    if queue is not None:
        qsub_cmd = "{} -q {}".format(qsub_cmd, queue)
    if stdout is not None:
        qsub_cmd = '{} -o {}'.format(qsub_cmd, stdout)
    if stderr is not None:
        qsub_cmd = '{} -e {}'.format(qsub_cmd, stderr)
    full_cmd = 'echo "{}" | {}'.format(mycmd, qsub_cmd)
    try:
        job_id = call(full_cmd)
        wait_for_job(job_id, sleep)
    finally:
        if job_state(job_id) not in 'CM':
            call('qdel {}'.format(job_id))
    if outfiles is not None:
        wait_for_files(outfiles, sleep, waittime)
    if not os.path.exists(sentinel):
        raise Exception('qsub sentinel not found: "{}"'.format(sentinel))
    with open(sentinel) as handle:
        code = handle.read().strip()
        if code != '0':
            raise Exception('qsub job "{}" exited with code "{}"'.format(full_cmd, code))


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
