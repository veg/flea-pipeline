#!/usr/bin/env python

"""
Runs the complete pipeline on a set of fasta files.

Input is a file with the following format:

<file_1> <id_1>
<file_2> <id_2>
...
<file_n> <id_n>


Usage:
  env_pipeline [options] <file>
  env_pipeline -h | --help

Options:
  --percent-identity=FLOAT  Percent identity for clustering [default: 0.99]
  --discard-lb=INT          Lower bound for selecting clusters [default: 5]
  --subsample-ub=INT        Upper bound for selecting clusters [default: 30]
  -v --verbose              Print progress.
  -h --help                 Show this screen.

"""

# TODO:
# - final alignment and trees
# - hardcoded paths
# - run jobs on cluster
# - update interface
# - for some reason clustering gets redone unnecessarily
# - order of inputs in tasks seems wrong sometimes
# - suppress job output
# - logging
# - update to latest versions of dependencies
# - replace unnecessary dependencies
# - update README

from subprocess import check_call
from subprocess import Popen, PIPE, STDOUT
import shlex
from os import path
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

#FIX THIS PLS
import os

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

from translate import translate
from backtranslate import backtranslate
from align_to_refs import align_to_refs, usearch_global
from DNAcons import dnacons
from correct_shifts import correct_shifts_fasta
from perfectORFs import perfect_file


def call(cmd_str, *args, stdout=None):
    """Call a command in the shell. Captures STDOUT and STDERR.

    If *args are given, they are passed as strings to STDIN.

    Raises an exception if return code != 0.

    """
    if stdout is None:
        stdout = PIPE
    if args:
        stdin = PIPE
    else:
        stdin = None
    process = Popen(cmd_str, stdin=stdin, stdout=stdout, stderr=PIPE, shell=True)
    if args:
        script_input = "".join("{}\n".format(i) for i in args)
        script_input = script_input.encode()
    else:
        script_input = None
    stdout_str, stderr_str = process.communicate(input=script_input)
    if process.returncode != 0:
        raise Exception("Failed to run '{}'\n{}{}Non-zero exit status {}".format(
                cmd_str, stdout_str, stderr_str, process.returncode))


def run_hyphy_script(script_file, *args, hyphy=None):
    if hyphy is None:
        hyphy = "HYPHYMP"
    cmd = '{} {}'.format(hyphy, script_file)
    call(cmd, *args)


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


def strlist(arg):
    if isinstance(arg, str):
        return [arg]
    return arg


def maybe(fun):
    @wraps(fun)
    def newf(infiles, outfiles):
        if any(os.stat(f).st_size == 0 for f in strlist(infiles)):
            for f in strlist(outfiles):
                touch(f)
        else:
            fun(infiles, outfiles)
    return newf


Timepoint = namedtuple('Timepoint', ['file', 'id', 'date'])


args = docopt(__doc__)

percent_identity = args["--percent-identity"]
discard_lb = int(args["--discard-lb"])
subsample_ub = int(args["--subsample-ub"])
is_verbose = args["--verbose"]

# FIXME: for now, filenames cannot have spaces. Make this more
# robust. For instance use tab-seperated values.

# FIXME: HYPHY does not like dots in sequence names.

with open(args['<file>'], newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    timepoints = list(Timepoint(f, i, d) for f, i, d in reader)

seq_ids = {t.file: t.id for t in timepoints}

start_files = list(t.file for t in timepoints)

data_dir = path.dirname(path.abspath(timepoints[0].file))
script_dir = os.path.split(__file__)[0]
# # make hyphy input and output directory
hyphy_data_dir = os.path.join(data_dir, "hyphy_data")
hyphy_input_dir = os.path.join(hyphy_data_dir, "input")
hyphy_results_dir = os.path.join(hyphy_data_dir, "results")

def _mkdir(d):
    try:
        os.mkdir(d)
    except FileExistsError:
        pass

_mkdir(hyphy_data_dir)
_mkdir(hyphy_input_dir)
_mkdir(hyphy_results_dir)

def hyphy_input(s):
    return os.path.join(hyphy_input_dir, s)


def hyphy_results(s):
    return os.path.join(hyphy_results_dir, s)


contaminant_db = '/home/keren/projects/env_pipeline/References/ContaminantRef.udb'
lanl_db = '/home/keren/projects/env_pipeline/References/LANLsubtypeRef.udb'


# @transform(start_files, suffix(".fastq"), '.clean.fastq')
# def clean_fastq(infile, outfile):
#     with open(outfile, 'w') as handle:
#         call('seqtk seq -l 0 {}'.format(infile), stdout=handle)


@transform(start_files, suffix(".fastq"), '.filtered.fasta')
def filter_fastq(infile, outfile):
    outfile = outfile[:-len('.fasta')]  # prinseq appends the extension
    call('prinseq -fastq {infile} -out_format 1 -out_good {outfile}'
         ' -min_len 2200 -max_len 2900 -min_qual_mean 25'
         ' -seq_id "{seq_id}_" -seq_id_mappings'.format(infile=infile,
                                                        outfile=outfile,
                                                        seq_id=seq_ids[infile]))


@transform(filter_fastq, suffix('.fasta'),
           ['.uncontaminated.fasta', '.contaminated.fasta'])
def filter_contaminants(infile, outfiles):
    uncontam, contam = outfiles
    call('usearch -usearch_global {infile} -db {db}'
         ' -id 0.98 -notmatched {uncontam} -matched {contam}'
         ' -strand both'.format(infile=infile, db=contaminant_db,
                                uncontam=uncontam, contam=contam))

@transform(filter_contaminants, suffix('.uncontaminated.fasta'),
           '.uncontaminated.matches.fasta')
def filter_lanl(infiles, outfile):
    uncontam, _ = infiles
    call('usearch -usearch_global {uncontam} -db {db} -id 0.8'
         ' -fastapairs {outfile} -strand both'.format(uncontam=uncontam,
                                                     db=lanl_db, outfile=outfile))


@transform(filter_lanl, suffix('.fasta'), '.shift-corrected.fasta')
def shift_correction(infile, outfile):
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile)
    n_dropped = n_seqs - n_fixed
    percent = 100 * n_dropped / n_seqs
    sumfile = '.'.join([outfile, 'summary'])
    with open(sumfile, 'w') as handle:
        handle.write('discarded {}/{} ({:.2f}%) sequences\n'.format(n_dropped, n_seqs, percent))


@transform(shift_correction, suffix('.fasta'), '.sorted.fasta')
def sort_by_length(infile, outfile):
    call('usearch -sortbylength {} -output {}'.format(infile, outfile))


@mkdir(sort_by_length, suffix('.fasta'), '.clusters')
@subdivide(sort_by_length, formatter(), '{path[0]}/{basename[0]}.clusters/cluster_*.raw.fasta')
def cluster(infile, outfiles, *args):
    for f in outfiles:
        os.unlink(f)
    outdir = '{}.clusters'.format(infile[:-len('.fasta')])
    outpattern = os.path.join(outdir, 'cluster_')
    call('usearch -cluster_smallmem {infile} -id {identity}'
         ' -clusters {outpattern}'.format(infile=infile, identity=percent_identity, outpattern=outpattern))

    # TODO: put this in own function
    r = re.compile(r'^cluster_[0-9]+$')
    for f in list(f for f in os.listdir(outdir) if r.match(f)):
        oldfile = os.path.join(outdir, f)
        newfile = ''.join([oldfile, '.raw.fasta'])
        os.rename(oldfile, newfile)


@transform(cluster, suffix('.fasta'), '.keep.fasta')
def select_cluster(infile, outfile):
    records = list(SeqIO.parse(infile, 'fasta'))
    if len(records) < discard_lb:
        # empty keep file; needed for pipeline sanity
        touch(outfile)
    if len(records) > subsample_ub:
        records = sample(records, subsample_ub)
    SeqIO.write(records, outfile, 'fasta')


@transform(select_cluster, suffix('.fasta'), '.aligned.fasta')
@maybe
def align_cluster(infile, outfile):
    with open(outfile, 'w') as handle:
        call('mafft {}'.format(infile), stdout=handle)


@transform(align_cluster, suffix('.fasta'), '.consensus.fasta')
@maybe
def cluster_consensus(infile, outfile):
    dnacons(infile, outfile=outfile, ungap=True)


@collate(cluster_consensus, formatter(), '{subdir[0][0]}.consensus.fasta')
def cat_clusters(infiles, outfile):
    cat_files(infiles, outfile)


@transform(cat_clusters, suffix('.fasta'), '.sorted.fasta')
def sort_consensus(infile, outfile):
    call('usearch -sortbylength {} -output {}'.format(infile, outfile))


@transform(sort_consensus, suffix('.fasta'), '.uniques.fasta')
def unique_consensus(infile, outfile):
    call('usearch -cluster_smallmem {} -id 1 -centroids {}'.format(infile, outfile))


@transform(unique_consensus, suffix('.fasta'), '.perfect.fasta')
def perfect_orfs(infile, outfile):
    perfect_file(infile, outfile, min_len=750, table=1, verbose=False)


@transform(perfect_orfs, suffix('.fasta'), '.udb')
def make_perfect_db(infile, outfile):
    call("usearch -makeudb_usearch {} -output {}".format(infile, outfile))


@collate([perfect_orfs, make_perfect_db, shift_correction],
         formatter(r'(?P<NAME>.+).filtered'),
         '{NAME[0]}.copynumber.fasta')
def add_copynumber(infiles, outfile):
    rawfile, perfectfile, dbfile = infiles
    pairs = usearch_global(rawfile, dbfile)
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


@merge(add_copynumber, "all_perfect_orfs.fasta")
def cat_all_perfect(infiles, outfile):
    cat_files(infiles, outfile)


@transform(cat_all_perfect, suffix('.fasta'), '.translated')
def translate_perfect(infile, outfile):
    translate(infile, outfile)


@transform(translate_perfect, formatter(), hyphy_input('merged.prot'))
def codon_align_perfect(infile, outfile):
    with open(outfile, 'w') as handle:
        call("mafft --auto {}".format(infile), stdout=handle)


@merge([cat_all_perfect, codon_align_perfect],
       hyphy_input('merged.fas'))
def backtranslate_alignment(infiles, outfile):
    aligned, perfect = infiles
    backtranslate(aligned, perfect, outfile)


@transform(backtranslate_alignment, formatter(), hyphy_input('merged.dates'))
def write_dates(infile, outfile):
    # NOTE: we assume the first part of the record id is the timestamp
    # id, followed by an underscore.
    records = list(SeqIO.parse(infile, "fasta"))
    id_to_date = {t.id : t.date for t in timepoints}
    with open(outfile, "w") as handle:
        outdict = {r.id: id_to_date[r.id.split("_")[0]] for r in records}
        json.dump(outdict, handle, separators=(",\n", ":"))


@transform(backtranslate_alignment, formatter(), hyphy_input('earlyCons.seq'))
def write_mrca(infile, outfile):
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(timepoints, key=strptime)
    oldest_records_filename = '.'.join([infile, "oldest_{}".format(oldest_timepoint.id)])
    mrca(infile, oldest_records_filename, outfile, oldest_timepoint.id)


@transform(backtranslate_alignment, formatter(), hyphy_input('merged.prot.parts'))
def write_hxb2_coords(infile, outfile):
    hxb2_script = os.path.join(script_dir, 'HXB2partsSplitter.bf')
    run_hyphy_script(hxb2_script, infile, outfile)


@merge([write_dates, write_hxb2_coords, backtranslate_alignment],
         [hyphy_input('mrca.seq')] + list(hyphy_results(f)
                                          for f in ('rates_pheno.tsv',
                                                    'trees.json',
                                                    'sequences.json')))
def evo_history(infiles, outfile):
    evolutionary_history_script = os.path.join(script_dir,
                                               "obtainEvolutionaryHistory.bf")
    run_hyphy_script(evolutionary_history_script, hyphy_data_dir)


@merge([write_dates, backtranslate_alignment, write_mrca],
       hyphy_results('frequences.json'))
def aa_freqs(infile, outfile):
    aa_freqs_script = os.path.join(script_dir,"aminoAcidFrequencies.bf")
    run_hyphy_script(aa_freqs_script, hyphy_data_dir)


@merge([write_dates, evo_history, backtranslate_alignment],
       hyphy_results('rates.json'))
def run_fubar(infile, outfile):
    fubar_script = os.path.join(script_dir, "runFUBAR.bf")
    run_hyphy_script(fubar_script, hyphy_data_dir)


pipeline_run([run_fubar, aa_freqs, evo_history], verbose=5, touch_files_only=True)

# # # run alignment in each timestep
# # timestep_aligned_files = []
# # for t in timepoints:
# #     # TODO: fix these names
# #     infile = "".join([t.file, ".pbformatfixed.fastq.good.fasta.matches.fasta.unshifted.fasta"])
# #     outfile = "".join([t.file, "_ALIGNED.fasta"])
# #     align_to_refs(infile, all_orfs_file, backtranslated_file,
# #                   db_file, outfile)
# #     timestep_aligned_files.append(outfile)

# # # concatenate all aligned timesteps
# # final_alignment_file = os.path.join(data_dir, "allTimepoints.ALIGNED.fasta")
# # cat_files(timestep_aligned_files, final_alignment_file)


# #FIX THIS HORRIBLE NONSENSE PLS
# # os.system("trees ALIGNED")
