#!/usr/bin/env python
"""
Runs the complete pipeline on a set of fasta files.

Input is a file with the following space-seperated informatoin on each
line for each time point:

<file> <id> <date>

Usage:
  env_pipeline [options] <file>
  env_pipeline -h | --help

"""

# TODO:
# - final alignment and trees
# - hardcoded paths
# - run jobs on cluster
# - logging
# - update to latest versions of dependencies
# - replace unnecessary dependencies
# - update README

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
from ruffus import cmdline

from translate import translate
from backtranslate import backtranslate
from align_to_refs import align_to_refs, usearch_global
from DNAcons import dnacons
from correct_shifts import correct_shifts_fasta
from perfectORFs import perfect_file
from util import call, hyphy_call, flatten, cat_files, touch, strlist


parser = cmdline.get_argparse(description='Run complete env pipeline')
parser.add_argument('file')
parser.add_argument('--id', default='0.99', type=str,
                    help='Percent identity for clustering.')
parser.add_argument('--min-cluster', default=5, type=int,
                    help='Minimum cluster size; smaller clusters are discarded.')
parser.add_argument('--max-cluster', default=30, type=int,
                    help='Maximum cluster size; larger clusters are subsampled.')

options = parser.parse_args()
if options.min_cluster < 1:
    raise Exception('Bad minimum cluster size: {}'. format(options.min_cluster))
if options.max_cluster < 1:
    raise Exception('Bad maximum cluster size: {}'. format(options.max_cluster))

# standard python logger which can be synchronised across concurrent
# Ruffus tasks
logger, logger_mutex = cmdline.setup_logging(__name__,
                                             options.log_file,
                                             options.verbose)


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


Timepoint = namedtuple('Timepoint', ['file', 'id', 'date'])


# FIXME: for now, filenames cannot have spaces. Make this more
# robust. For instance use tab-seperated values.

# FIXME: HYPHY does not like dots in sequence names.

with open(options.file, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    timepoints = list(Timepoint(f, i, d) for f, i, d in reader)

seq_ids = {t.file: t.id for t in timepoints}
start_files = list(t.file for t in timepoints)

data_dir = path.dirname(path.abspath(timepoints[0].file))
script_dir = os.path.split(__file__)[0]
hyphy_data_dir = os.path.join(data_dir, "hyphy_data")
hyphy_input_dir = os.path.join(hyphy_data_dir, "input")
hyphy_results_dir = os.path.join(hyphy_data_dir, "results")


def hyphy_input(s):
    return os.path.join(hyphy_input_dir, s)


def hyphy_results(s):
    return os.path.join(hyphy_results_dir, s)


contaminant_db = '/home/keren/projects/env_pipeline/References/ContaminantRef.udb'
lanl_db = '/home/keren/projects/env_pipeline/References/LANLsubtypeRef.udb'


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
         ' -clusters {outpattern}'.format(infile=infile, identity=options.id,
                                          outpattern=outpattern))

    # TODO: put this in own function
    r = re.compile(r'^cluster_[0-9]+$')
    for f in list(f for f in os.listdir(outdir) if r.match(f)):
        oldfile = os.path.join(outdir, f)
        newfile = ''.join([oldfile, '.raw.fasta'])
        os.rename(oldfile, newfile)


@transform(cluster, suffix('.fasta'), '.keep.fasta')
def select_cluster(infile, outfile):
    records = list(SeqIO.parse(infile, 'fasta'))
    if len(records) < options.min_cluster:
        # empty keep file; needed for pipeline sanity
        touch(outfile)
    if len(records) > options.max_cluster:
        records = sample(records, options.max_cluster)
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
    assert(perfectfile[-len('.perfect.fasta'):] == '.perfect.fasta')
    assert(dbfile[-len('.udb'):] == '.udb')
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


@mkdir(hyphy_input_dir)
@mkdir(hyphy_results_dir)
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
    perfect, aligned = infiles
    assert(os.path.basename(perfect) == 'all_perfect_orfs.fasta')
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
    hyphy_call(hxb2_script, infile, outfile)


@merge([write_dates, write_hxb2_coords, backtranslate_alignment],
         [hyphy_input('mrca.seq')] + list(hyphy_results(f)
                                          for f in ('rates_pheno.tsv',
                                                    'trees.json',
                                                    'sequences.json')))
def evo_history(infiles, outfile):
    evolutionary_history_script = os.path.join(script_dir,
                                               "obtainEvolutionaryHistory.bf")
    hyphy_call(evolutionary_history_script, hyphy_data_dir)


@merge([write_dates, backtranslate_alignment, write_mrca],
       hyphy_results('frequences.json'))
def aa_freqs(infile, outfile):
    aa_freqs_script = os.path.join(script_dir,"aminoAcidFrequencies.bf")
    hyphy_call(aa_freqs_script, hyphy_data_dir)


@merge([write_dates, evo_history, backtranslate_alignment],
       hyphy_results('rates.json'))
def run_fubar(infile, outfile):
    fubar_script = os.path.join(script_dir, "runFUBAR.bf")
    hyphy_call(fubar_script, hyphy_data_dir)


cmdline.run(options)

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
