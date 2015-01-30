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
from align_to_refs import align_to_refs, usearch_global
from DNAcons import dnacons
from correct_shifts import correct_shifts_fasta, write_correction_result
from perfectORFs import perfect_file
from util import call, qsub, hyphy_call, cat_files, touch, strlist, traverse
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

seq_ids = {t.file: t.id for t in timepoints}
start_files = list(t.file for t in timepoints)

for f in start_files:
    if not os.path.exists(f):
        raise Exception('file does not exist: "{}"'.format(f))


if len(seq_ids) != len(timepoints):
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
        raise Exception('Sequnce ratios do not match. Expected {}:{}. Got: {}:{}. Inputs: "{}" Output: "{}"'.format(a_exp, b_exp, a_n, b_n, inputs, output))


def must_work(maybe=False, seq_ratio=None, seq_ids=False, pattern=None):
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
        return wrapped
    return wrap


def check_suffix(name, suffix):
    assert(name.endswith(suffix))


def remove_suffix(name, suffix):
    check_suffix(name, suffix)
    return name[:-len(suffix)]


def check_basename(name, bn):
    assert(os.path.basename(name) == bn)


@originate('this_run.config')
def write_config(outfile):
    with open(outfile, 'w') as handle:
        config.write(handle)


@jobs_limit(n_local_jobs, local_job_limiter)
@follows(write_config)
@transform(start_files, suffix(".fastq"), '.filtered.fasta')
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
            min_qual_mean=min_qual_mean, seq_id=seq_ids[infile]))


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(filter_fastq, suffix('.fasta'),
           ['.uncontaminated.fasta', '.contaminated.fasta'])
def filter_contaminants(infile, outfiles):
    uncontam, contam = outfiles
    call('{usearch} -usearch_global {infile} -db {db}'
         ' -id {id} -notmatched {uncontam} -matched {contam}'
         ' -strand both'.format(usearch=config['Paths']['usearch'],
                                infile=infile, db=config['Paths']['contaminants_db'],
                                id=config['Parameters']['contaminant_identity'],
                                uncontam=uncontam, contam=contam))


def db_search_pairs(infile, outfile, dbfile, identity, nums_only=False):
    if nums_only:
        cmd = ("{usearch} -usearch_global {infile} -db {db} -id {id}"
               " -userout {outfile} -userfields query+target -strand both")
    else:
        cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
               ' -fastapairs {outfile} -strand both')
    call(cmd.format(usearch=config['Paths']['usearch'],
                    infile=infile, db=dbfile, id=identity, outfile=outfile))


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(filter_contaminants, suffix('.uncontaminated.fasta'),
           '.uncontaminated.matches.fasta')
@must_work()
def filter_uncontaminated(infiles, outfile):
    uncontam, _ = infiles
    db_search_pairs(uncontam, outfile, config['Paths']['reference_db'],
                    config['Parameters']['reference_identity'])


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(filter_uncontaminated, suffix('.fasta'), '.shift-corrected.fasta')
@must_work(seq_ratio=(2, 1))
def shift_correction(infile, outfile):
    correct_shifts_fasta(infile, outfile, keep=True)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(shift_correction, suffix('.fasta'), '.sorted.fasta')
@must_work(seq_ids=True)
def sort_by_length(infile, outfile):
    call('{usearch} -sortbylength {infile} -output {outfile}'.format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))


cluster_flag = "cluster_task_completed.flag"


def cluster_uptodate(infile, outfiles):
    """A temporary hack until the issue with @subdivide and @split
    gets resolved.

    """
    if not os.path.exists(cluster_flag):
        return True, "{} is missing".format(cluster_flag)
    # ensure each sorted file has a cluster directory
    sorted_files = glob('*.sorted.fasta')
    sorted_files = list(f for f in sorted_files
                        if not f.endswith('ref_pairs.shift-corrected.sorted.fasta'))
    for f in sorted_files:
        d = '{}.clusters'.format(remove_suffix(f, '.fasta'))
        if not os.path.exists(d):
            return True, "at least one cluster directory is missing: {}".format(d)
    pattern = '{}.clusters/*.raw.fasta'
    # ensure each cluster directory contains at least one cluster
    for f in sorted_files:
        search = pattern.format(remove_suffix(f, '.fasta'))
        raw_clusters = glob(search)
        if not raw_clusters:
            return True, "at least one cluster directory is empty: {}".format(d)
    # ensure all timestamps are more recent than sorted file timestamps.
    for f in sorted_files:
        sorted_time = os.path.getmtime(f)
        search = pattern.format(remove_suffix(f, '.fasta'))
        raw_clusters = glob(search)
        for r in raw_clusters:
            if os.path.getmtime(r) <= sorted_time:
                return True, "at least one cluster is old: {}".format(r)
    return False, "{} exists and clusters look up-to-date".format(cluster_flag)


@check_if_uptodate(cluster_uptodate)
@posttask(touch_file(cluster_flag))
@jobs_limit(n_local_jobs, local_job_limiter)
@mkdir(sort_by_length, suffix('.fasta'), '.clusters')
@subdivide(sort_by_length, formatter(), '{path[0]}/{basename[0]}.clusters/cluster_*.raw.fasta')
def cluster(infile, outfiles):
    for f in outfiles:
        os.unlink(f)
    outdir = '{}.clusters'.format(infile[:-len('.fasta')])
    outpattern = os.path.join(outdir, 'cluster_')
    call('{usearch} -cluster_smallmem {infile} -id {id}'
         ' -clusters {outpattern}'.format(
            usearch=config['Paths']['usearch'],
            infile=infile, id=config['Parameters']['cluster_identity'],
            outpattern=outpattern))

    # TODO: put this in own function
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


def rm_std():
    """Get rid of annoying STDIN files' after qsub runs"""
    for f in glob('STDIN.*'):
        os.unlink(f)


def maybe_qsub(cmd, sentinel):
    if use_cluster:
        if os.path.exists(sentinel):
            os.unlink(sentinel)
        qsub(cmd, sentinel, walltime=int(config['Misc']['align_clusters_walltime']))
    else:
        call(cmd)


@posttask(rm_std)
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(select_clusters, suffix('.fasta'), '.aligned.fasta')
@must_work(maybe=True)
def align_clusters(infile, outfile):
    stderr = '{}.stderr'.format(outfile)
    sentinel = '{}.complete'.format(outfile)
    cmd = 'mafft --quiet {} > {} 2>{}'.format(infile, outfile, stderr)
    maybe_qsub(cmd, sentinel)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(align_clusters, suffix('.fasta'), '.consensus.fasta')
@must_work(maybe=True)
def cluster_consensus(infile, outfile):
    dnacons(infile, outfile=outfile, ungap=True)


@jobs_limit(n_local_jobs, local_job_limiter)
@collate(cluster_consensus, formatter(), '{subdir[0][0]}.consensus.fasta')
@must_work(seq_ids=True)
def cat_clusters(infiles, outfile):
    cat_files(infiles, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(cat_clusters, suffix('.fasta'), '.ref_pairs.fasta')
@must_work()
def consensus_db_search(infile, outfile):
    db_search_pairs(infile, outfile, config['Paths']['reference_db'],
                    config['Parameters']['reference_identity'])


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(consensus_db_search, suffix('.fasta'), '.shift-corrected.fasta')
@must_work()
def consensus_shift_correction(infile, outfile):
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile, keep=False)
    sumfile = '{}.summary'.format(outfile)
    write_correction_result(n_seqs, n_fixed, sumfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(consensus_shift_correction, suffix('.fasta'), '.sorted.fasta')
@must_work(seq_ids=True)
def sort_consensus(infile, outfile):
    call('{usearch} -sortbylength {infile} -output {outfile}'.format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(sort_consensus, suffix('.fasta'), '.uniques.fasta')
@must_work()
def unique_consensus(infile, outfile):
    call('{usearch} -cluster_smallmem {infile} -id 1 -centroids {outfile}'.format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(unique_consensus, suffix('.fasta'), '.perfect.fasta')
@must_work()
def perfect_orfs(infile, outfile):
    perfect_file(infile, outfile, min_len=int(config['Parameters']['min_orf_length']),
                 table=1, verbose=False)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(perfect_orfs, suffix('.fasta'), '.udb')
@must_work()
def make_individual_dbs(infile, outfile):
    call("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))


@jobs_limit(n_local_jobs, local_job_limiter)
@collate([perfect_orfs, make_individual_dbs, shift_correction],
         formatter(r'(?P<NAME>.+).filtered'),
         '{NAME[0]}.copynumber.fasta')
@must_work(seq_ratio=(1, 1), pattern='*perfect*fasta')
def add_copynumber(infiles, outfile):
    rawfile, perfectfile, dbfile = infiles
    check_suffix(perfectfile, '.perfect.fasta')
    check_suffix(dbfile, '.udb')
    identity = config['Parameters']['raw_to_consensus_identity']
    pairs = usearch_global(rawfile, dbfile, identity, usearch=config['Paths']['usearch'])
    consensus_counts = defaultdict(lambda: 1)
    for raw_id, ref_id in pairs:
        consensus_counts[str(ref_id).upper()] += 1
    records = list(SeqIO.parse(perfectfile, 'fasta'))
    def record_generator():
        for r in records:
            str_id = str(r.id).upper()
            r.id = "{}_{}".format(str_id, consensus_counts[str_id])
            r.name = ""
            r.description = ""
            yield r
    SeqIO.write(record_generator(), outfile, 'fasta')


@mkdir(hyphy_input_dir)
@mkdir(hyphy_results_dir)
@jobs_limit(n_local_jobs, local_job_limiter)
@merge(add_copynumber, "all_perfect_orfs.fasta")
@must_work(seq_ids=True)
def cat_all_perfect(infiles, outfile):
    cat_files(infiles, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(cat_all_perfect, suffix('.fasta'), '.udb')
@must_work()
def make_full_db(infile, outfile):
    call("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=config['Paths']['usearch'], infile=infile, outfile=outfile))


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(cat_all_perfect, suffix('.fasta'), '.translated')
@must_work(seq_ids=True)
def translate_perfect(infile, outfile):
    translate(infile, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(translate_perfect, formatter(), hyphy_input('merged.prot'))
@must_work(seq_ids=True)
def codon_align_perfect(infile, outfile):
    with open(outfile, 'w') as handle:
        call("mafft --auto {}".format(infile), stdout=handle)


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
    dnacons(recordfile, id_str="mrca", ungap=False, outfile=outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(backtranslate_alignment, formatter(), hyphy_input('earlyCons.seq'))
@must_work()
def compute_mrca(infile, outfile):
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(timepoints, key=strptime)
    oldest_records_filename = '.'.join([infile, "oldest_{}".format(oldest_timepoint.id)])
    mrca(infile, oldest_records_filename, outfile, oldest_timepoint.id)


@active_if(config.getboolean('Tasks', 'hyphy'))
@jobs_limit(n_local_jobs, local_job_limiter)
@transform(backtranslate_alignment, formatter(), hyphy_input('merged.prot.parts'))
@must_work()
def compute_hxb2_coords(infile, outfile):
    hxb2_script = hyphy_script('HXB2partsSplitter.bf')
    hyphy_call(hxb2_script, infile, outfile)


@active_if(config.getboolean('Tasks', 'hyphy'))
@jobs_limit(n_local_jobs, local_job_limiter)
@merge([write_dates, compute_hxb2_coords, backtranslate_alignment],
         [hyphy_input('mrca.seq')] + list(hyphy_results(f)
                                          for f in ('rates_pheno.tsv',
                                                    'trees.json',
                                                    'sequences.json')))
@must_work()
def evo_history(infiles, outfile):
    hyphy_call(hyphy_script("obtainEvolutionaryHistory.bf"), hyphy_data_dir)


@active_if(config.getboolean('Tasks', 'hyphy'))
@jobs_limit(n_local_jobs, local_job_limiter)
@merge([write_dates, backtranslate_alignment, compute_mrca],
       hyphy_results('frequencies.json'))
@must_work()
def aa_freqs(infile, outfile):
    hyphy_call(hyphy_script("aminoAcidFrequencies.bf"), hyphy_data_dir)


@active_if(config.getboolean('Tasks', 'hyphy'))
@jobs_limit(n_local_jobs, local_job_limiter)
@merge([write_dates, evo_history, backtranslate_alignment],
       hyphy_results('rates.json'))
@must_work()
def run_fubar(infile, outfile):
    hyphy_call(hyphy_script('runFUBAR.bf'), hyphy_data_dir)


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_local_jobs, local_job_limiter)
@transform(shift_correction,
           suffix('.shift-corrected.fasta'),
           add_inputs([ make_full_db]),
           '.shift-corrected.raw_consensus_pairs.txt')
@must_work()
def full_timestep_pairs(infiles, outfile):
    infile, (dbfile,) = infiles
    check_suffix(infile, '.shift-corrected.fasta')
    check_suffix(dbfile, '.udb')
    identity = config['Parameters']['reference_identity']
    db_search_pairs(infile, outfile, dbfile, identity, nums_only=True)


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_local_jobs, local_job_limiter)
@subdivide(full_timestep_pairs,
           formatter('.*/(?P<NAME>.+).raw_consensus_pairs'),
           add_inputs([cat_all_perfect]),
           '{NAME[0]}.raw_consensus_pairs.combined.*.unaligned.fasta',
           '{NAME[0]}')
@must_work()
def combine_pairs(infiles, outfiles, basename):
    for f in outfiles:
        os.unlink(f)
    infile, (perfectfile,) = infiles
    readsfile = '{}.fasta'.format(basename)
    check_basename(perfectfile, 'all_perfect_orfs.fasta')
    with open(infile) as handle:
        pairs = list(line.strip().split("\t") for line in handle.readlines())
    match_dict = defaultdict(list)
    for read, ref in pairs:
        match_dict[ref].append(read)
    references_records = list(SeqIO.parse(perfectfile, "fasta"))
    read_records = list(SeqIO.parse(readsfile, "fasta"))
    references_dict = {r.id : r for r in references_records}
    read_dict = {r.id : r for r in read_records}

    for ref_id, read_ids in match_dict.items():
        outfile = '{}.raw_consensus_pairs.combined.{}.unaligned.fasta'.format(basename, ref_id)
        records = [references_dict[ref_id]]
        records.extend(list(read_dict[i] for i in read_ids))
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
    maybe_qsub(cmd, sentinel)


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(codon_align, suffix('.bam'), '.fasta')
@must_work()
def convert_bam_to_fasta(infile, outfile):
    cmd = "bam2msa.py {} {}".format(infile, outfile)
    sentinel = '{}.complete'.format(outfile)
    maybe_qsub(cmd, sentinel)


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_local_jobs, local_job_limiter)
@transform(convert_bam_to_fasta,
           suffix('.fasta'),
           add_inputs([backtranslate_alignment]),
           '.msa.fasta')
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
