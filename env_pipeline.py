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
# - check everywhere usearch is used; summarize # seqs lost
# - are identities re-used correctly?
# - decorator to ensure output not empty
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

from translate import translate
from backtranslate import backtranslate
from align_to_refs import align_to_refs, usearch_global
from DNAcons import dnacons
from correct_shifts import correct_shifts_fasta, write_correction_result
from perfectORFs import perfect_file
from util import call, qsub, hyphy_call, flatten, cat_files, touch, strlist


parser = cmdline.get_argparse(description='Run complete env pipeline',
                              ignored_args=['jobs'])
parser.add_argument('file')
parser.add_argument('--config', type=str,
                    help='Configuration file.')

options = parser.parse_args()

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
max_n_jobs = max(n_local_jobs, n_remote_jobs)

options.jobs = max_n_jobs


def hyphy_script(name):
    return os.path.join(hyphy_script_dir, name)


def hyphy_input(s):
    return os.path.join(hyphy_input_dir, s)


def hyphy_results(s):
    return os.path.join(hyphy_results_dir, s)


def maybe(fun):
    """Modify a task to create empty output files if any input file is
    empty."""
    @wraps(fun)
    def newf(infiles, outfiles):
        if any(os.stat(f).st_size == 0 for f in strlist(infiles)):
            for f in strlist(outfiles):
                touch(f)
        else:
            fun(infiles, outfiles)
    return newf


def check_suffix(name, suffix):
    assert(name[-len(suffix):] == suffix)


def check_basename(name, bn):
    assert(os.path.basename(name) == bn)


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(start_files, suffix(".fastq"), '.filtered.fasta')
def filter_fastq(infile, outfile):
    outfile = outfile[:-len('.fasta')]  # prinseq appends the extension
    call('prinseq -fastq {infile} -out_format 1 -out_good {outfile}'
         ' -min_len 2200 -max_len 2900 -min_qual_mean 25'
         ' -seq_id "{seq_id}_" -seq_id_mappings'.format(infile=infile,
                                                        outfile=outfile,
                                                        seq_id=seq_ids[infile]))


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(filter_fastq, suffix('.fasta'),
           ['.uncontaminated.fasta', '.contaminated.fasta'])
def filter_contaminants(infile, outfiles):
    uncontam, contam = outfiles
    call('usearch -usearch_global {infile} -db {db}'
         ' -id {id} -notmatched {uncontam} -matched {contam}'
         ' -strand both'.format(infile=infile, db=config['Paths']['contaminants_db'],
                                id=config['Parameters']['contaminant_identity'],
                                uncontam=uncontam, contam=contam))


def db_search_pairs(infile, outfile, dbfile, identity):
    call('usearch -usearch_global {infile} -db {db} -id {id}'
         ' -fastapairs {outfile} -strand both'.format(
            infile=infile, db=dbfile, id=identity, outfile=outfile))


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(filter_contaminants, suffix('.uncontaminated.fasta'),
           '.uncontaminated.matches.fasta')
def filter_uncontaminated(infiles, outfile):
    uncontam, _ = infiles
    db_search_pairs(uncontam, outfile, config['Paths']['reference_db'],
                    config['Parameters']['reference_identity'])


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(filter_uncontaminated, suffix('.fasta'), '.shift-corrected.fasta')
def shift_correction(infile, outfile):
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile, keep=True)
    sumfile = '{}.summary'.format(outfile)
    write_correction_result(n_seqs, n_fixed, sumfile)


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(shift_correction, suffix('.fasta'), '.sorted.fasta')
def sort_by_length(infile, outfile):
    call('usearch -sortbylength {} -output {}'.format(infile, outfile))


@jobs_limit(n_local_jobs, 'local_jobs')
@mkdir(sort_by_length, suffix('.fasta'), '.clusters')
@subdivide(sort_by_length, formatter(), '{path[0]}/{basename[0]}.clusters/cluster_*.raw.fasta')
def cluster(infile, outfiles, *args):
    for f in outfiles:
        os.unlink(f)
    outdir = '{}.clusters'.format(infile[:-len('.fasta')])
    outpattern = os.path.join(outdir, 'cluster_')
    call('usearch -cluster_smallmem {infile} -id {id}'
         ' -clusters {outpattern}'.format(infile=infile,
                                          id=config['Parameters']['cluster_identity'],
                                          outpattern=outpattern))

    # TODO: put this in own function
    r = re.compile(r'^cluster_[0-9]+$')
    for f in list(f for f in os.listdir(outdir) if r.match(f)):
        oldfile = os.path.join(outdir, f)
        newfile = ''.join([oldfile, '.raw.fasta'])
        os.rename(oldfile, newfile)


@jobs_limit(n_local_jobs, 'local_jobs')
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


@jobs_limit(n_remote_jobs, 'remote_jobs')
@transform(select_clusters, suffix('.fasta'), '.aligned.fasta')
@maybe
def align_clusters(infile, outfile):
    sentinel = '{}.job.complete'.format(outfile)
    stderr = '{}.job.stderr'.format(outfile)
    if os.path.exists(sentinel):
        os.unlink(sentinel)
    # --quiet option is needed; otherwise it fails
    qsub('mafft --quiet {} > {} 2>{}'.format(infile, outfile, stderr), sentinel)
    statinfo = os.stat(outfile)
    if statinfo.st_size == 0:
        raise Exception('mafft produced empty output file: "{}"'.format(outfile))


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(align_clusters, suffix('.fasta'), '.consensus.fasta')
@maybe
def cluster_consensus(infile, outfile):
    dnacons(infile, outfile=outfile, ungap=True)


@jobs_limit(n_local_jobs, 'local_jobs')
@collate(cluster_consensus, formatter(), '{subdir[0][0]}.consensus.fasta')
def cat_clusters(infiles, outfile):
    cat_files(infiles, outfile)
    statinfo = os.stat(outfile)
    if statinfo.st_size == 0:
        raise Exception('A timepoint has no consensus sequences:'
                        ' "{}"'.format(outfile))


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(cat_clusters, suffix('.fasta'), '.ref_pairs.fasta')
def consensus_db_search(infile, outfile):
    db_search_pairs(infile, outfile, config['Paths']['reference_db'],
                    config['Parameters']['reference_identity'])


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(consensus_db_search, suffix('.fasta'), '.shift-corrected.fasta')
def consensus_shift_correction(infile, outfile):
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile, keep=False)
    sumfile = '{}.summary'.format(outfile)
    write_correction_result(n_seqs, n_fixed, sumfile)


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(consensus_shift_correction, suffix('.fasta'), '.sorted.fasta')
def sort_consensus(infile, outfile):
    call('usearch -sortbylength {} -output {}'.format(infile, outfile))


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(sort_consensus, suffix('.fasta'), '.uniques.fasta')
def unique_consensus(infile, outfile):
    call('usearch -cluster_smallmem {} -id 1 -centroids {}'.format(infile, outfile))


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(unique_consensus, suffix('.fasta'), '.perfect.fasta')
def perfect_orfs(infile, outfile):
    perfect_file(infile, outfile, min_len=750, table=1, verbose=False)


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(perfect_orfs, suffix('.fasta'), '.udb')
def make_individual_dbs(infile, outfile):
    call("usearch -makeudb_usearch {} -output {}".format(infile, outfile))


@jobs_limit(n_local_jobs, 'local_jobs')
@collate([perfect_orfs, make_individual_dbs, shift_correction],
         formatter(r'(?P<NAME>.+).filtered'),
         '{NAME[0]}.copynumber.fasta')
def add_copynumber(infiles, outfile):
    rawfile, perfectfile, dbfile = infiles
    check_suffix(perfectfile, '.perfect.fasta')
    check_suffix(dbfile, '.udb')
    identity = config['Parameters']['reference_identity']
    pairs = usearch_global(rawfile, dbfile, identity)
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
@jobs_limit(n_local_jobs, 'local_jobs')
@merge(add_copynumber, "all_perfect_orfs.fasta")
def cat_all_perfect(infiles, outfile):
    cat_files(infiles, outfile)


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(cat_all_perfect, suffix('.fasta'), '.udb')
def make_full_db(infile, outfile):
    call("usearch -makeudb_usearch {} -output {}".format(infile, outfile))


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(cat_all_perfect, suffix('.fasta'), '.translated')
def translate_perfect(infile, outfile):
    translate(infile, outfile)


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(translate_perfect, formatter(), hyphy_input('merged.prot'))
def codon_align_perfect(infile, outfile):
    with open(outfile, 'w') as handle:
        call("mafft --auto {}".format(infile), stdout=handle)


@jobs_limit(n_local_jobs, 'local_jobs')
@merge([cat_all_perfect, codon_align_perfect],
       hyphy_input('merged.fas'))
def backtranslate_alignment(infiles, outfile):
    perfect, aligned = infiles
    check_basename(perfect, 'all_perfect_orfs.fasta')
    backtranslate(aligned, perfect, outfile)


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(backtranslate_alignment, formatter(), hyphy_input('merged.dates'))
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


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(backtranslate_alignment, formatter(), hyphy_input('earlyCons.seq'))
def compute_mrca(infile, outfile):
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(timepoints, key=strptime)
    oldest_records_filename = '.'.join([infile, "oldest_{}".format(oldest_timepoint.id)])
    mrca(infile, oldest_records_filename, outfile, oldest_timepoint.id)


@jobs_limit(n_local_jobs, 'local_jobs')
@transform(backtranslate_alignment, formatter(), hyphy_input('merged.prot.parts'))
def compute_hxb2_coords(infile, outfile):
    hxb2_script = hyphy_script('HXB2partsSplitter.bf')
    hyphy_call(hxb2_script, infile, outfile)


@jobs_limit(n_local_jobs, 'local_jobs')
@merge([write_dates, compute_hxb2_coords, backtranslate_alignment],
         [hyphy_input('mrca.seq')] + list(hyphy_results(f)
                                          for f in ('rates_pheno.tsv',
                                                    'trees.json',
                                                    'sequences.json')))
def evo_history(infiles, outfile):
    hyphy_call(hyphy_script("obtainEvolutionaryHistory.bf"), hyphy_data_dir)


@jobs_limit(n_local_jobs, 'local_jobs')
@merge([write_dates, backtranslate_alignment, compute_mrca],
       hyphy_results('frequences.json'))
def aa_freqs(infile, outfile):
    hyphy_call(hyphy_script("aminoAcidFrequencies.bf"), hyphy_data_dir)


@jobs_limit(n_local_jobs, 'local_jobs')
@merge([write_dates, evo_history, backtranslate_alignment],
       hyphy_results('rates.json'))
def run_fubar(infile, outfile):
    hyphy_call(hyphy_script('runFUBAR.bf'), hyphy_data_dir)


@active_if(config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_local_jobs, 'local_jobs')
@transform(shift_correction,
           suffix('.shift-corrected.fasta'),
           add_inputs([cat_all_perfect, backtranslate_alignment, make_full_db]),
           '.shift-corrected.aligned.fasta')
def align_full_timestep(infiles, outfile):
    raise Exception('Update the pipeline code to run HYPHY scripts on the output'
                    ' of this task before running it.')
    shift_corrected, (perfect, perfect_aligned, dbfile) = infiles

    check_suffix(shift_corrected, '.shift-corrected.fasta')
    check_basename(perfect, 'all_perfect_orfs.fasta')
    check_suffix(dbfile, '.udb')
    check_basename(perfect_aligned, 'merged.fas')
    identity = config['Parameters']['reference_identity']
    align_to_refs(shift_corrected, perfect, perfect_aligned, dbfile,
                  outfile, identity)


@jobs_limit(n_local_jobs, 'local_jobs')
@active_if(config.getboolean('Tasks', 'align_full'))
@merge(align_full_timestep, 'all_timepoints.aligned.fasta')
def merge_all_timepoints(infiles, outfile):
    cat_files(infiles, outfile)


@jobs_limit(n_local_jobs, 'local_jobs')
@active_if(config.getboolean('Tasks', 'align_full') and
           config.getboolean('Tasks', 'generate_trees'))
@transform([align_full_timestep, merge_all_timepoints],
           suffix('.fasta'), '.tre')
def infer_trees(infile, outfile):
    with open(infile, 'rb') as in_handle:
        with open(outfile, 'w') as out_handle:
            call('FastTreeMP -nt -gtr -nosupport', stdin=in_handle, stdout=out_handle)


if __name__ == '__main__':
    cmdline.run(options)
