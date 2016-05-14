import os
import re
from random import sample
import tempfile
from functools import partial
from collections import defaultdict
import itertools
import warnings
import shutil

import numpy as np
from ruffus import Pipeline, suffix, formatter, add_inputs
from Bio import SeqIO

from flea_pipeline.util import maybe_qsub, cat_files, touch, strlist, traverse
from flea_pipeline.util import new_record_seq_str, insert_gaps, update_record_seq
from flea_pipeline.util import must_work, must_produce, report_wrapper
from flea_pipeline.util import local_job_limiter, remote_job_limiter
from flea_pipeline.util import check_suffix, check_basename, n_jobs
from flea_pipeline.util import read_single_record
from flea_pipeline.util import partition
from flea_pipeline.util import grouper
from flea_pipeline.util import run_regexp
from flea_pipeline.util import translate_helper
from flea_pipeline.util import remove_suffix

import flea_pipeline.pipeline_globals as globals_


pipeline_dir = os.path.join(globals_.data_dir, "quality")


def min_len():
    min_len_fraction = globals_.config.getfloat('Parameters', 'min_sequence_length')
    ref_file = globals_.config.get('Parameters', 'reference_sequence')
    # multiple by 3 because reference sequence is translated.
    return 3 * int(min_len_fraction * len(read_single_record(ref_file, 'fasta').seq))


@must_work()
@report_wrapper
def filter_fastq(infile, outfile):
    qmax = globals_.config.get('Parameters', 'qmax')
    max_err_rate = globals_.config.get('Parameters', 'max_err_rate')
    cmd = ('{usearch} -fastq_filter {infile} -fastq_maxee_rate {max_err_rate}'
         ' -threads 1 -fastq_qmax {qmax} -fastq_minlen {min_len} -fastqout {outfile}'
         ' -relabel "{seq_id}_"'.format(
            usearch=globals_.config.get('Paths', 'usearch'),
            infile=infile, outfile=outfile, qmax=qmax, min_len=min_len(),
            max_err_rate=max_err_rate, seq_id=globals_.key_to_label[infile]))
    return maybe_qsub(cmd, infile, outfiles=outfile, name='filter-fastq')


@must_work()
@report_wrapper
def trim_heads_and_tails(infile, outfile):
    python = globals_.config.get('Paths', 'python')
    script = os.path.join(globals_.script_dir, 'trim.py')
    cmd = ("{python} {script} --fastq {infile} {outfile}")
    cmd = cmd.format(python=python, script=script, infile=infile, outfile=outfile)
    return maybe_qsub(cmd, infile, outfile, name='trim-heads-and-tails')


@must_work()
@report_wrapper
def filter_runs(infile, outfile):
    runlen = globals_.config.getint('Parameters', 'run_length')
    cregexp = run_regexp(runlen)
    records = SeqIO.parse(infile, 'fastq')
    result = (r for r in records if cregexp.search(str(r.seq)) is None)
    SeqIO.write(result, outfile, 'fastq')


@must_work()
@report_wrapper
def filter_contaminants(infile, outfile):
    uncontam = outfile
    contam = outfile.replace('uncontam', 'contam')
    db = globals_.config.get('Parameters', 'contaminants_db')
    id_ = globals_.config.get('Parameters', 'contaminant_identity')
    outfiles = [uncontam, contam]
    cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
           ' -notmatchedfq {uncontam} -matchedfq {contam}'
           ' -strand both'.format(usearch=globals_.config.get('Paths', 'usearch'),
                                  infile=infile, db=db, id=id_,
                                  uncontam=uncontam, contam=contam))
    return maybe_qsub(cmd, infile, outfiles=outfiles, name='filter-contaminants')


@must_work()
@report_wrapper
def filter_db(infile, outfile):
    """Keep sequences that match the database.

    Also produces a 'calns' file with strand information (+ or -) and
    compressed pairwise alignment information for each sequence and
    its match.
    """
    dbfile = globals_.config.get('Parameters', 'reference_db')
    identity = globals_.config.get('Parameters', 'reference_identity')
    max_accepts = globals_.config.get('Parameters', 'max_accepts')
    max_rejects = globals_.config.get('Parameters', 'max_rejects')
    useroutfile = '{}.userout'.format(remove_suffix(outfile, '.fastq'))
    outfiles = [outfile, useroutfile]

    cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
           ' -matchedfq {outfile} -userout {userout}'
           ' -userfields qstrand+tstrand+caln'
           ' -top_hit_only -strand both'
           ' -maxaccepts {max_accepts} -maxrejects {max_rejects}')
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, db=dbfile, id=identity,
                     outfile=outfile, userout=useroutfile,
                     max_accepts=max_accepts, max_rejects=max_rejects)
    return maybe_qsub(cmd, infile, outfiles=outfile, name="filter-db")


@must_work()
@report_wrapper
def trim_terminal_gaps(infile, outfile):
    """Get reverse complement of each sequence if necessary. Then trim
    terminal insertions relative to the matched database sequence.

    """
    userfile = '{}.userout'.format(remove_suffix(infile, '.fastq'))
    infiles = [infile, userfile]
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "trim_terminal_gaps.py"),
        'infile': infile,
        'userfile': userfile,
        'outfile': outfile,
        }
    cmd = ("{python} {script} {infile} {userfile} {outfile}").format(**kwargs)
    return maybe_qsub(cmd, infiles, outfile, name="trim-terminal-gaps")


@must_work()
@report_wrapper
def filter_length(infile, outfile):
    cutoff = min_len()
    records = SeqIO.parse(infile, 'fastq')
    result = (r for r in records if len(r.seq) >= cutoff)
    SeqIO.write(result, outfile, 'fastq')


def make_quality_pipeline(name=None):
    if name is None:
        name = "quality_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    filter_fastq_task = pipeline.transform(filter_fastq,
                                           input=None,
                                           filter=formatter(),
                                           output=os.path.join(pipeline_dir, '{basename[0]}.qfilter.fastq'))
    filter_fastq_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    trim_task = pipeline.transform(trim_heads_and_tails,
                                    input=filter_fastq_task,
                                    filter=suffix('.fastq'),
                                    output='.trimmed.fastq')
    trim_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    filter_runs_task = pipeline.transform(filter_runs,
                                          input=trim_task,
                                          filter=suffix('.fastq'),
                                          output='.noruns.fastq')
    filter_runs_task.jobs_limit(n_local_jobs, local_job_limiter)

    filter_contaminants_task = pipeline.transform(filter_contaminants,
                                                  input=filter_runs_task,
                                                  filter=suffix('.fastq'),
                                                  output='.uncontam.fastq')
    filter_contaminants_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    filter_db_task = pipeline.transform(filter_db,
                                        input=filter_contaminants_task,
                                        filter=suffix('.fastq'),
                                        output='.dbfiltered.fastq')
    filter_db_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    trim_terminal_gaps_task = pipeline.transform(trim_terminal_gaps,
                                            input=filter_db_task,
                                            filter=suffix('.fastq'),
                                            output='.gapstripped.fastq')
    trim_terminal_gaps_task.jobs_limit(n_local_jobs, local_job_limiter)

    filter_length_task = pipeline.transform(filter_length,
                                            input=trim_terminal_gaps_task,
                                            filter=suffix('.fastq'),
                                            output='.lenfilter.fastq')
    filter_length_task.jobs_limit(n_local_jobs, local_job_limiter)

    pipeline.set_head_tasks([filter_fastq_task])
    pipeline.set_tail_tasks([filter_length_task])

    for task in pipeline.head_tasks:
        task.mkdir(pipeline_dir)

    return pipeline