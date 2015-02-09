#!/usr/bin/env python
"""
Runs the complete pipeline on a set of fasta files.

Input is a file with the following space-seperated informatoin on each
line for each time point:

<file> <id> <date>

Usage:
  env_pipeline.py [options] <file>
  env_pipeline.py -h | --help

"""

# TODO:
# - run as many tasks as possible on the cluster
# - order of input files not guaranteed; make more robust
# - check everywhere usearch is used; summarize # seqs lost
# - are identities re-used correctly?
# - rename files; flea and test subdirectories
# - logging
# - update to latest versions of dependencies
# - replace unnecessary dependencies
# - update README

# FIXME: for now, filenames cannot have spaces. Make this more
# robust. For instance use tab-seperated values.

# FIXME: HYPHY does not like dots in sequence names.


import shlex
import os
import csv
import json
from collections import namedtuple, defaultdict
from datetime import datetime
import sys
import uuid
import shutil
import re
from random import sample
from functools import wraps
from configparser import ConfigParser, ExtendedInterpolation
from glob import glob
from fnmatch import fnmatch
import tempfile

from docopt import docopt

from Bio import SeqIO

from ruffus import pipeline_run
from ruffus import suffix
from ruffus import regex
from ruffus import transform
from ruffus import subdivide
from ruffus import collate
from ruffus import mkdir
from ruffus import formatter
from ruffus import merge
from ruffus import cmdline
from ruffus import add_inputs
from ruffus import active_if
from ruffus import jobs_limit
from ruffus import check_if_uptodate
from ruffus import touch_file
from ruffus import posttask
from ruffus import originate
from ruffus import follows


from translate import translate
from backtranslate import backtranslate
from DNAcons import dnacons
from correct_shifts import correct_shifts_fasta, write_correction_result
from perfectORFs import perfect_file
from util import call, qsub, cat_files, touch, strlist, traverse
from util import new_record_seq_str, insert_gaps


parser = cmdline.get_argparse(description='Run complete env pipeline',
                              ignored_args=['jobs'])
parser.add_argument('file')
parser.add_argument('--config', type=str,
                    help='Configuration file.')

options = parser.parse_args()

if not options.verbose:
    options.verbose = 1

# standard python logger which can be synchronised across concurrent
# Ruffus tasks
logger, logger_mutex = cmdline.setup_logging(__name__,
                                             options.log_file,
                                             options.verbose)


Timepoint = namedtuple('Timepoint', ['file', 'id', 'date'])

with open(options.file, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    timepoints = list(Timepoint(f, i, d) for f, i, d in reader)

timepoint_ids = {t.file: t.id for t in timepoints}
start_files = list(t.file for t in timepoints)

for f in start_files:
    if not os.path.exists(f):
        raise Exception('file does not exist: "{}"'.format(f))


if len(set(t.id for t in timepoints)) != len(timepoints):
    raise Exception('non-unique sequence ids')

# useful directories
data_dir = os.path.dirname(os.path.abspath(timepoints[0].file))
script_dir = os.path.abspath(os.path.split(__file__)[0])
hyphy_script_dir = os.path.join(script_dir, 'hyphy_scripts')
hyphy_data_dir = os.path.join(data_dir, "hyphy_data")
hyphy_input_dir = os.path.join(hyphy_data_dir, "input")
hyphy_results_dir = os.path.join(hyphy_data_dir, "results")

config = ConfigParser(interpolation=ExtendedInterpolation())

if options.config is None:
    configfile = os.path.join(script_dir, 'env_pipeline.config')
else:
    configfile = options.config
config.read(configfile)


n_local_jobs = int(config['Jobs']['local_jobs'])
n_remote_jobs = int(config['Jobs']['remote_jobs'])
use_cluster = config.getboolean('Jobs', 'use_cluster')

if n_local_jobs < 1:
    n_local_jobs = 1

if n_remote_jobs < 1 and use_cluster:
    raise Exception('Bad parameters; use_cluster="{}"'
                    ' but remote_jobs="{}"'.format(use_cluster, n_remote_jobs))

if not use_cluster:
    n_remote_jobs = n_local_jobs

options.jobs = max(n_local_jobs, n_remote_jobs)

local_job_limiter = 'local_jobs'
remote_job_limiter = 'remote_jobs'


def hyphy_script(name):
    return os.path.join(hyphy_script_dir, name)


def hyphy_input(s):
    return os.path.join(hyphy_input_dir, s)


def hyphy_results(s):
    return os.path.join(hyphy_results_dir, s)


def ensure_not_empty(files):
    for f in strlist(files):
        if not os.path.exists(f):
            raise Exception('Expected output file does not'
                            ' exist: "{}"'.format(f))
        if os.stat(f).st_size == 0:
            raise Exception('Empty file: "{}"'.format(f))


def check_seq_ids(inputs, output):
    a_ids = set(r.id for a in inputs for r in SeqIO.parse(a, 'fasta'))
    b_ids = set(r.id for r in SeqIO.parse(output, 'fasta'))
    if a_ids != b_ids:
        raise Exception('IDs from "{}" do not match "{}"'.format(inputs, output))


def check_seq_ratio(inputs, output, expected):
    a_exp, b_exp = expected
    assert a_exp == int(a_exp)
    assert b_exp == int(b_exp)
    a_n = sum(1 for a in inputs for r in SeqIO.parse(a, 'fasta'))
    b_n = sum(1 for r in SeqIO.parse(output, 'fasta'))
    if a_n * b_exp != b_n * a_exp:
        raise Exception('Sequnce ratios do not match.'
                        ' Expected {}:{}. Got: {}:{}.'
                        ' Inputs: "{}" Output: "{}"'.format(
                a_exp, b_exp, a_n, b_n, inputs, output))


def check_illegal_chars(f, chars):
    for r in SeqIO.parse(f, 'fasta'):
        found = set(chars).intersection(set(str(r.seq)))
        if found:
            raise Exception('Illegal characters "{}" found in sequence "{}"'
                            ' of file "{}"'.format(found, r.id, f))


def must_work(maybe=False, seq_ratio=None, seq_ids=False, illegal_chars=None, pattern=None):
    """Fail if any output is empty.

    maybe: touch output and return if any input is empty
    seq_ratio: (int, int): ensure numbers of sequences match
    seq_ids: ensure all ids present in input files are in output file
    pattern: pattern for determing input files to consider

    """
    if pattern is None:
        pattern = '*'
    def wrap(fun):
        if options.touch_files_only or options.just_print:
            return fun
        @wraps(fun)  # necessary because ruffus uses function name internally
        def wrapped(infiles, outfiles, *args, **kwargs):
            if maybe:
                if any(os.stat(f).st_size == 0 for f in traverse(strlist(infiles))):
                    for f in traverse(strlist(outfiles)):
                        touch(f)
                    return
            fun(infiles, outfiles, *args, **kwargs)
            infiles = list(traverse(strlist(infiles)))
            infiles = list(f for f in infiles if fnmatch(f, pattern))
            outfiles = strlist(outfiles)
            ensure_not_empty(outfiles)
            if seq_ids:
                assert len(outfiles) == 1
                check_seq_ids(infiles, outfiles[0])
            if seq_ratio is not None:
                assert len(outfiles) == 1
                check_seq_ratio(infiles, outfiles[0], seq_ratio)
            if illegal_chars:
                for f in outfiles:
                    check_illegal_chars(f, illegal_chars)
        return wrapped
    return wrap


def must_produce(n):
    """Check that at least `n` files match the pathname glob"""
    def wrap(fun):
        @wraps(fun)
        def wrapped(infiles, outfiles, pathname, *args, **kwargs):
            fun(infiles, outfiles, pathname, *args, **kwargs)
            n_produced = len(glob(pathname))
            if n_produced < n:
                raise Exception('Task was supposed to produce at least {}'
                                ' outputs, but it produced {}'.format(n, n_produced))
        return wrapped
    return wrap


def maybe_qsub(cmd, sentinel, **kwargs):
    if use_cluster:
        if os.path.exists(sentinel):
            os.unlink(sentinel)
        qsub(cmd, sentinel, walltime=int(config['Misc']['align_clusters_walltime']),
             queue=config['Jobs']['queue'], nodes=config['Jobs']['nodes'],
             ppn=config['Jobs']['ppn'], **kwargs)
    else:
        # TODO: handle stdout and stderr kwargs
        call(cmd)


def mafft(infile, outfile):
    stderr = '{}.stderr'.format(outfile)
    sentinel = '{}.complete'.format(outfile)
    cmd = 'mafft --quiet {} > {} 2>{}'.format(infile, outfile, stderr)
    maybe_qsub(cmd, sentinel, outfiles=outfile, stdout='/dev/null', stderr='/dev/null')


def check_suffix(name, suffix):
    assert(name.endswith(suffix))


def remove_suffix(name, suffix):
    check_suffix(name, suffix)
    return name[:-len(suffix)]


def check_basename(name, bn):
    assert(os.path.basename(name) == bn)


def param_uptodate(infile, outfile):
    """Re-run if any parameters changed."""
    # TODO: do this check on a task basis
    if not os.path.exists(outfile):
        return True, "{} is missing".format(outfile)
    param_config = ConfigParser(interpolation=ExtendedInterpolation())
    param_config.read(outfile)
    if param_config['Parameters'] != config['Parameters']:
        return True, 'parameters not the same'
    return False, 'parameters are the same'


@check_if_uptodate(param_uptodate)
@originate('parameters.config')
def write_config(outfile):
    with open(outfile, 'w') as handle:
        config.write(handle)


@jobs_limit(n_local_jobs, local_job_limiter)
@follows(write_config)
@transform(start_files, suffix(".fastq"), '.qfilter.fasta')
@must_work()  # my decorators must go before ruffus ones
def filter_fastq(infile, outfile):
    outfile = outfile[:-len('.fasta')]  # prinseq appends the extension
    min_len = config['Parameters']['min_sequence_length']
    max_len = config['Parameters']['max_sequence_length']
    min_qual_mean = config['Parameters']['min_qual_mean']
    call('{prinseq} -fastq {infile} -out_format 1 -out_good {outfile}'
         ' -min_len {min_len} -max_len {max_len} -min_qual_mean {min_qual_mean}'
         ' -seq_id "{seq_id}_" -seq_id_mappings'.format(
            prinseq=config['Paths']['prinseq'],
            infile=infile, outfile=outfile, min_len=min_len, max_len=max_len,
            min_qual_mean=min_qual_mean, seq_id=timepoint_ids[infile]))


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(filter_fastq, suffix('.fasta'),
           ['.uncontam.fasta', '.contam.fasta'])
def filter_contaminants(infile, outfiles):
    uncontam, contam = outfiles
    sentinel = '{}.complete'.format(contam)
    cmd = ('{usearch} -usearch_global {infile} -db {db}'
           ' -id {id} -notmatched {uncontam} -matched {contam}'
           ' -strand both'.format(usearch=config['Paths']['usearch'],
                                  infile=infile, db=config['Paths']['contaminants_db'],
                                  id=config['Parameters']['contaminant_identity'],
                                  uncontam=uncontam, contam=contam))
    maybe_qsub(cmd, sentinel, outfiles=outfiles, name='filter-contaminants')


def usearch_global_pairs(infile, outfile, dbfile, identity, nums_only=False, name=None):
    if nums_only:
        cmd = ("{usearch} -usearch_global {infile} -db {db} -id {id}"
               " -userout {outfile} -userfields query+target -strand both")
    else:
        cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
               ' -fastapairs {outfile} -strand both')
    sentinel = '{}.complete'.format(outfile)
    if name is None:
        name = 'usearch-global-pairs'
    cmd = cmd.format(usearch=config['Paths']['usearch'],
                     infile=infile, db=dbfile, id=identity, outfile=outfile)
    maybe_qsub(cmd, sentinel, outfiles=outfile, name=name)


def usearch_global_get_pairs(infile, dbfile, identity, name=None):
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, 'pairs.txt')
        usearch_global_pairs(infile, outfile, dbfile, identity, nums_only=True, name=name)
        with open(outfile) as f:
            return list(line.strip().split("\t") for line in f.readlines())


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(filter_contaminants, suffix('.uncontam.fasta'),
           '.uncontam.rfilter.fasta')
@must_work()
def filter_uncontaminated(infiles, outfile):
    uncontam, _ = infiles
    usearch_global_pairs(uncontam, outfile, config['Paths']['reference_db'],
                         config['Parameters']['reference_identity'],
                         name="filter-uncontaminated")


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(filter_uncontaminated, suffix('.fasta'), '.shifted.fasta')
@must_work(seq_ratio=(2, 1))
def shift_correction(infile, outfile):
    correct_shifts_fasta(infile, outfile, keep=True)


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(shift_correction, suffix('.fasta'), '.sorted.fasta')
@must_work(seq_ids=True)
def sort_by_length(infile, outfile):
    sentinel = "{}.complete".format(outfile)
    cmd = ('{usearch} -sortbylength {infile} -output {outfile}'.format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))
    maybe_qsub(cmd, sentinel, outfiles=outfile, name="sort-by-length")


@jobs_limit(n_remote_jobs, remote_job_limiter)
@mkdir(sort_by_length, suffix('.fasta'), '.clusters')
@subdivide(sort_by_length, formatter(),
           '{path[0]}/{basename[0]}.clusters/cluster_*.raw.fasta',
           '{path[0]}/{basename[0]}.clusters/cluster_*.raw.fasta')
@must_produce(n=int(config['Parameters']['min_n_clusters']))
def cluster(infile, outfiles, pathname):
    for f in outfiles:
        os.unlink(f)
    outdir = '{}.clusters'.format(infile[:-len('.fasta')])
    outpattern = os.path.join(outdir, 'cluster_')
    cmd = ('{usearch} -cluster_smallmem {infile} -id {id}'
           ' -clusters {outpattern}'.format(
            usearch=config['Paths']['usearch'],
            infile=infile, id=config['Parameters']['cluster_identity'],
            outpattern=outpattern))
    sentinel = "cluster.complete"
    maybe_qsub(cmd, sentinel, name="cluster")
    r = re.compile(r'^cluster_[0-9]+$')
    for f in list(f for f in os.listdir(outdir) if r.match(f)):
        oldfile = os.path.join(outdir, f)
        newfile = ''.join([oldfile, '.raw.fasta'])
        os.rename(oldfile, newfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(cluster, suffix('.fasta'), '.keep.fasta')
def select_clusters(infile, outfile):
    records = list(SeqIO.parse(infile, 'fasta'))
    minsize = int(config['Parameters']['min_cluster_size'])
    maxsize = int(config['Parameters']['max_cluster_size'])
    if len(records) < minsize:
        touch(outfile)  # empty file; needed for pipeline sanity
        return
    if len(records) > maxsize:
        records = sample(records, maxsize)
    SeqIO.write(records, outfile, 'fasta')


def hyphy_call(script_file, name, args):
    if args:
        in_str = "".join("{}\n".format(i) for i in args)
    else:
        in_str = ''
    infile = '{}.stdin'.format(name)
    sentinel = '{}.complete'.format(name)
    with open(infile, 'w') as handle:
        handle.write(in_str)
    cmd = '{} {} < {}'.format(config['Paths']['hyphy'], script_file, infile)
    maybe_qsub(cmd, sentinel, name=name)


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(select_clusters, suffix('.fasta'), '.aligned.fasta')
@must_work(maybe=True)
def align_clusters(infile, outfile):
    mafft(infile, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(align_clusters, suffix('.fasta'), '.cons.fasta')
@must_work(maybe=True, illegal_chars='X-')
def cluster_consensus(infile, outfile):
    dnacons(infile, outfile, ungap=True)


@jobs_limit(n_local_jobs, local_job_limiter)
@collate(cluster_consensus, formatter(), '{subdir[0][0]}.cons.fasta')
@must_work(seq_ids=True)
def cat_clusters(infiles, outfile):
    cat_files(infiles, outfile)


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(cat_clusters, suffix('.fasta'), '.pairs.fasta')
@must_work()
def consensus_db_search(infile, outfile):
    usearch_global_pairs(infile, outfile, config['Paths']['reference_db'],
                         config['Parameters']['reference_identity'],
                         name='consensus-db-search')


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(consensus_db_search, suffix('.fasta'), '.shifted.fasta')
@must_work()
def consensus_shift_correction(infile, outfile):
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile, keep=False)
    sumfile = '{}.summary'.format(outfile)
    write_correction_result(n_seqs, n_fixed, sumfile)


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(consensus_shift_correction, suffix('.fasta'), '.sorted.fasta')
@must_work(seq_ids=True)
def sort_consensus(infile, outfile):
    cmd = ('{usearch} -sortbylength {infile} -output {outfile}'.format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))
    sentinel = '{}.complete'.format(outfile)
    maybe_qsub(cmd, sentinel, outfiles=outfile, name='sort-consensus')


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(sort_consensus, suffix('.fasta'), '.uniques.fasta')
@must_work()
def unique_consensus(infile, outfile):
    sentinel = '{}.complete'.format(outfile)
    cmd = ('{usearch} -cluster_smallmem {infile} -id 1 -centroids {outfile}'.format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))
    maybe_qsub(cmd, sentinel, outfiles=outfile, name='unique_consensus')


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(unique_consensus, suffix('.fasta'), '.perfect.fasta')
@must_work()
def perfect_orfs(infile, outfile):
    perfect_file(infile, outfile, min_len=int(config['Parameters']['min_orf_length']),
                 table=1, verbose=False)


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(perfect_orfs, suffix('.fasta'), '.udb')
@must_work()
def make_individual_dbs(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))
    sentinel = '{}.complete'.format(outfile)
    maybe_qsub(cmd, sentinel, outfiles=outfile, name='make-individual-dbs')


@jobs_limit(n_remote_jobs, remote_job_limiter)
@collate([perfect_orfs, make_individual_dbs, shift_correction],
         formatter(r'(?P<NAME>.+).qfilter'),
         '{NAME[0]}.copynumber.fasta',
         '{NAME[0]}')
@must_work(seq_ratio=(1, 1), pattern='*perfect*fasta')
def add_copynumber(infiles, outfile, basename):
    rawfile, perfectfile, dbfile = infiles
    check_suffix(perfectfile, '.perfect.fasta')
    check_suffix(dbfile, '.udb')
    identity = config['Parameters']['raw_to_consensus_identity']
    pairfile = '{}.copynumber.pairs'.format(basename)
    usearch_global_pairs(rawfile, pairfile, dbfile, identity,
                         nums_only=True, name='add-copynumber')
    with open(pairfile) as f:
            pairs = list(line.strip().split("\t") for line in f.readlines())
    consensus_counts = defaultdict(lambda: 1)
    for raw_id, ref_id in pairs:
        consensus_counts[str(ref_id).upper()] += 1
    def new_record(r):
        r = r[:]
        _id = str(r.id).upper()
        r.id = "{}_{}".format(_id, consensus_counts[_id])
        r.name = ''
        r.description = ''
        return r
    SeqIO.write((new_record(r) for r in SeqIO.parse(perfectfile, 'fasta')),
                outfile, 'fasta')


@jobs_limit(n_local_jobs, local_job_limiter)
@merge(add_copynumber, "all_perfect_orfs.fasta")
@must_work(seq_ids=True)
def cat_all_perfect(infiles, outfile):
    cat_files(infiles, outfile)


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(cat_all_perfect, suffix('.fasta'), '.udb')
@must_work()
def make_full_db(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))
    sentinel = '{}.complete'.format(outfile)
    maybe_qsub(cmd, sentinel, outfiles=outfile, name='make_full_db')


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(cat_all_perfect, suffix('.fasta'), '.translated.fasta')
@must_work(seq_ids=True)
def translate_perfect(infile, outfile):
    translate(infile, outfile)


@mkdir(hyphy_input_dir)
@mkdir(hyphy_results_dir)
@jobs_limit(n_local_jobs, local_job_limiter)
@transform(translate_perfect, formatter(), hyphy_input('merged.prot'))
@must_work(seq_ids=True)
def codon_align_perfect(infile, outfile):
    mafft(infile, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@merge([cat_all_perfect, codon_align_perfect],
       hyphy_input('merged.fas'))
@must_work(seq_ids=True)
def backtranslate_alignment(infiles, outfile):
    perfect, aligned = infiles
    check_basename(perfect, 'all_perfect_orfs.fasta')
    backtranslate(aligned, perfect, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(backtranslate_alignment, formatter(), hyphy_input('merged.dates'))
@must_work()
def write_dates(infile, outfile):
    # NOTE: we assume the first part of the record id is the timestamp
    # id, followed by an underscore.
    records = list(SeqIO.parse(infile, "fasta"))
    id_to_date = {t.id : t.date for t in timepoints}
    with open(outfile, "w") as handle:
        outdict = {r.id: id_to_date[r.id.split("_")[0]] for r in records}
        json.dump(outdict, handle, separators=(",\n", ":"))


def mrca(infile, recordfile, outfile, oldest_id):
    """Writes records from `infile` to `recordfile` that have an ID
    corresponding to `oldest_id`. Then runs DNAcons, writing result to
    `outfile`.

    """
    len_id = len(oldest_id)
    records = SeqIO.parse(infile, "fasta")
    oldest_records = (r for r in records if r.id.startswith(oldest_id))
    SeqIO.write(oldest_records, recordfile, "fasta")
    dnacons(recordfile, outfile, id_str="mrca", ungap=False, ambiguous='N')


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(backtranslate_alignment, formatter(), hyphy_input('earlyCons.seq'))
@must_work(illegal_chars='X')
def compute_mrca(infile, outfile):
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(timepoints, key=strptime)
    oldest_records_filename = '.'.join([infile, "oldest_{}".format(oldest_timepoint.id)])
    mrca(infile, oldest_records_filename, outfile, oldest_timepoint.id)


@active_if(config.getboolean('Tasks', 'hyphy'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(backtranslate_alignment, formatter(), hyphy_input('merged.prot.parts'))
@must_work()
def compute_hxb2_coords(infile, outfile):
    hxb2_script = hyphy_script('HXB2partsSplitter.bf')
    hyphy_call(hxb2_script, 'hxb2_coords', [infile, outfile])


@active_if(config.getboolean('Tasks', 'hyphy'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@merge([write_dates, compute_hxb2_coords, backtranslate_alignment, compute_mrca],
         [hyphy_input('mrca.seq')] + list(hyphy_results(f)
                                          for f in ('rates_pheno.tsv',
                                                    'trees.json',
                                                    'sequences.json')))
@must_work()
def evo_history(infiles, outfile):
    hyphy_call(hyphy_script("obtainEvolutionaryHistory.bf"),
               'evo_history', [hyphy_data_dir])


@active_if(config.getboolean('Tasks', 'hyphy'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@merge([write_dates, backtranslate_alignment, compute_mrca],
       hyphy_results('frequencies.json'))
@must_work()
def aa_freqs(infile, outfile):
    hyphy_call(hyphy_script("aminoAcidFrequencies.bf"),
               'aa_freqs',
               [hyphy_data_dir])


@active_if(config.getboolean('Tasks', 'hyphy'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@merge([write_dates, evo_history, backtranslate_alignment],
       hyphy_results('rates.json'))
@must_work()
def run_fubar(infile, outfile):
    hyphy_call(hyphy_script('runFUBAR.bf'), 'fubar', [hyphy_data_dir])


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(shift_correction,
           suffix('.fasta'),
           add_inputs([ make_full_db]),
           '.pairs.txt')
@must_work()
def full_timestep_pairs(infiles, outfile):
    infile, (dbfile,) = infiles
    identity = config['Parameters']['raw_to_consensus_identity']
    usearch_global_pairs(infile, outfile, dbfile, identity, nums_only=True,
                         name='full-timestep-pairs')


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_local_jobs, local_job_limiter)
@mkdir(full_timestep_pairs, suffix('.pairs.txt'), '.alignments')
@subdivide(full_timestep_pairs,
           formatter('.*/(?P<NAME>.+).pairs.txt'),
           add_inputs([cat_all_perfect]),
           '{NAME[0]}.alignments/combined.*.unaligned.fasta',
           '{NAME[0]}')
@must_work()
def combine_pairs(infiles, outfiles, basename):
    for f in outfiles:
        os.unlink(f)
    infile, (perfectfile,) = infiles
    seqsfile = '{}.fasta'.format(basename)
    check_basename(perfectfile, 'all_perfect_orfs.fasta')
    with open(infile) as handle:
        pairs = list(line.strip().split("\t") for line in handle.readlines())
    match_dict = defaultdict(list)
    for seq, ref in pairs:
        match_dict[ref].append(seq)
    references_records = list(SeqIO.parse(perfectfile, "fasta"))
    seq_records = list(SeqIO.parse(seqsfile, "fasta"))
    references_dict = {r.id : r for r in references_records}
    seq_dict = {r.id : r for r in seq_records}

    for ref_id, seq_ids in match_dict.items():
        outfile = '{}.alignments/combined.{}.unaligned.fasta'.format(basename, ref_id)
        records = [references_dict[ref_id]]
        records.extend(list(seq_dict[i] for i in seq_ids))
        SeqIO.write(records, outfile, "fasta")


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(combine_pairs,
           suffix('.unaligned.fasta'),
           '.aligned.bam')
@must_work()
def codon_align(infile, outfile):
    cmd = "bealign.py -R {} {}".format(infile, outfile)
    sentinel = '{}.complete'.format(outfile)
    stdout = '{}.stdout'.format(outfile)
    stderr = '{}.stderr'.format(outfile)
    maybe_qsub(cmd, sentinel, outfiles=outfile, stdout=stdout, stderr=stderr)


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(codon_align, suffix('.bam'), '.fasta')
@must_work()
def convert_bam_to_fasta(infile, outfile):
    cmd = "bam2msa.py {} {}".format(infile, outfile)
    sentinel = '{}.complete'.format(outfile)
    stdout = '{}.stdout'.format(outfile)
    stderr = '{}.stderr'.format(outfile)
    maybe_qsub(cmd, sentinel, outfiles=outfile, stdout=stdout, stderr=stderr)


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_local_jobs, local_job_limiter)
@transform(convert_bam_to_fasta,
           suffix('.fasta'),
           add_inputs([backtranslate_alignment]),
           '.gapped.fasta')
@must_work()
def insert_gaps_task(infiles, outfile):
    infile, (backtranslated,) = infiles
    ref, *seqs = list(SeqIO.parse(infile, 'fasta'))
    ref_gapped = list(r for r in SeqIO.parse(backtranslated, 'fasta')
                      if r.id == ref.id)[0]
    seqs_gapped = (new_record_seq_str(r, insert_gaps(str(ref_gapped.seq),
                                                     str(r.seq),
                                                     '-', '-', skip=1))
                   for r in seqs)
    SeqIO.write(seqs_gapped, outfile, "fasta")


@jobs_limit(n_local_jobs, local_job_limiter)
@active_if(config.getboolean('Tasks', 'align_full'))
@merge(insert_gaps_task, 'all_timepoints.aligned.fasta')
@must_work(seq_ids=True)
def merge_all_timepoints(infiles, outfile):
    cat_files(infiles, outfile)


if __name__ == '__main__':
    checksum_level = int(config['Misc']['checksum_level'])
    cmdline.run(options, checksum_level=checksum_level)
