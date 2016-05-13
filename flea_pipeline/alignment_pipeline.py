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


pipeline_dir = os.path.join(globals_.data_dir, "alignment")


def mafft(infile, outfile):
    stderr = '{}.stderr'.format(outfile)
    cmd = 'mafft-fftns --ep 0.5 --quiet --preservecase {} > {} 2>{}'.format(infile, outfile, stderr)
    return maybe_qsub(cmd, infile, outfiles=outfile)


@must_work()
@report_wrapper
def translate_wrapper(infile, outfile):
    return translate_helper(infile, outfile, gapped=False, name='translate')


@must_work(seq_ids=True)
@report_wrapper
def mafft_wrapper_seq_ids(infile, outfile):
    mafft(infile, outfile)


@must_work(seq_ids=True)
@report_wrapper
def cat_all_hqcs(infiles, outfile):
    if len(infiles) != len(globals_.timepoints):
        raise Exception('Number of input files does not match number'
                        ' of timepoints')
    cat_files(infiles, outfile)


@must_work()
@report_wrapper
def copy_file(infile, outfile):
    shutil.copyfile(infile, outfile)


def pause_for_editing_alignment():
    if globals_.config.getboolean('Tasks', 'pause_after_codon_alignment'):
        infile = os.path.join(pipeline_dir, 'hqcs.translated.aligned.fasta')
        outfile = os.path.join(pipeline_dir, 'hqcs.translated.aligned.edited.fasta')
        input('Paused for manual editing of {}'
              '\nPress Enter to continue.'.format(outfile))
        # check that all sequences in `outfile` differ only in gap locations
        unedited = dict((r.id, r) for r in SeqIO.parse(infile, 'fasta'))
        edited = dict((r.id, r) for r in SeqIO.parse(outfile, 'fasta'))
        for k in edited.keys():
            edited_seq = list(c for c in str(edited[k].seq) if c != '-')
            unedited_seq = list(c for c in str(unedited[k].seq) if c != '-')
            if edited_seq != unedited_seq:
                raise Exception('Incorrect manual editing of "{}".'
                                ' "{}" differs after editing.'.format(outfile, k))


@must_work()
@report_wrapper
def backtranslate_alignment(infiles, outfile):
    hqcs, aligned_protein = infiles
    check_basename(hqcs, 'hqcs.fasta')
    check_suffix(aligned_protein, '.translated.aligned.edited.fasta')
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "backtranslate.py"),
        'protein': aligned_protein,
        'dna': hqcs,
        'outfile': outfile,
        }
    cmd = "{python} {script} {protein} {dna} {outfile}".format(**kwargs)
    stdout = os.path.join(globals_.qsub_dir, 'backtranslate.stdout')
    stderr = os.path.join(globals_.qsub_dir, 'backtranslate.stderr')
    return maybe_qsub(cmd, infiles, outfiles=outfile,
                      stdout=stdout, stderr=stderr,
                      name="backtranslate")


def make_alignment_pipeline(name=None):
    if name is None:
        name = "alignment_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    translate_hqcs_task = pipeline.transform(translate_wrapper,
                                             name='translate_hqcs',
                                             input=None,
                                             filter=formatter('.*/(?P<LABEL>.+).fasta'),
                                             output=os.path.join(pipeline_dir, '{LABEL[0]}.translated.fasta'))
    translate_hqcs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    align_hqcs_protein_task = pipeline.transform(mafft_wrapper_seq_ids,
                                                 name='align_hqcs_protein',
                                                 input=translate_hqcs_task,
                                                 filter=suffix('.fasta'),
                                                 output='.aligned.fasta')
    align_hqcs_protein_task.jobs_limit(n_local_jobs, local_job_limiter)

    copy_protein_alignment_task = pipeline.transform(copy_file,
                                                     name='copy_protein_alignment',
                                                     input=align_hqcs_protein_task,
                                                     filter=suffix('.fasta'),
                                                     output='.edited.fasta')
    copy_protein_alignment_task.jobs_limit(n_local_jobs, local_job_limiter)
    copy_protein_alignment_task.posttask(pause_for_editing_alignment)

    # expects all HQCSs as input
    backtranslate_alignment_task = pipeline.transform(backtranslate_alignment,
                                                      input=None,
                                                      add_inputs=add_inputs(copy_protein_alignment_task),
                                                      filter=formatter(),
                                                      output=os.path.join(pipeline_dir, 'hqcs.translated.aligned.edited.backtranslated.fasta'))
    backtranslate_alignment_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    pipeline.set_head_tasks([translate_hqcs_task, backtranslate_alignment_task])
    pipeline.set_tail_tasks([backtranslate_alignment_task, copy_protein_alignment_task])

    for task in pipeline.head_tasks:
        task.mkdir(pipeline_dir)

    return pipeline
