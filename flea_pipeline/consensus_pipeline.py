import os
import re
import random
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
from flea_pipeline.util import usearch_hqcs_ids

import flea_pipeline.pipeline_globals as globals_


pipeline_dir = os.path.join(globals_.data_dir, "consensus")


@must_work()
@report_wrapper
def cluster(infile, outfile):
    minsl = globals_.config.get('Parameters', 'min_length_ratio')
    cmd = ('{usearch} -cluster_fast {infile} -id {id}'
           ' -uc {outfile} -sort length'
           ' -maxaccepts {max_accepts} -maxrejects {max_rejects}'
           ' -top_hit_only -minsl {minsl}')
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, id=globals_.config.get('Parameters', 'cluster_identity'),
                     outfile=outfile,
                     max_accepts=globals_.config.get('Parameters', 'max_accepts'),
                     max_rejects=globals_.config.get('Parameters', 'max_rejects'),
                     minsl=minsl)
    return maybe_qsub(cmd, infile, outfile, name="cluster")


@must_work()
@report_wrapper
def fastq_clusters(infiles, outfiles, outdir):
    for f in outfiles:
        os.unlink(f)
    fastqfile, ucfile = infiles
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "cluster_fastq.py"),
        'ucfile': ucfile,
        'fastqfile': fastqfile,
        'outdir': outdir,
        'minsize': globals_.config.get('Parameters', 'min_cluster_size'),
        }
    cmd = ("{python} {script} --minsize {minsize}"
           " {ucfile} {fastqfile} {outdir}".format(**kwargs))
    return maybe_qsub(cmd, infiles, [], name="cluster-fastq")


def iter_sample(iterable, k):
    """Sample `k` items from an iterable.

    Adapted from http://stackoverflow.com/a/12581484

    """
    results = []
    for i, v in enumerate(iterable):
        r = random.randint(0, i)
        if r < k:
            if i < k:
                results.insert(r, v) # add first `n` items in random order
            else:
                results[r] = v # at a decreasing rate, replace random items
    if len(results) < k:
        raise ValueError("Sample larger than population.")
    return results


@must_work()
@report_wrapper
def sample_clusters(infile, outfile):
    n = sum(1 for r in SeqIO.parse(infile, 'fastq'))
    records = SeqIO.parse(infile, 'fastq')
    maxsize = int(globals_.config.get('Parameters', 'max_cluster_size'))
    if n > maxsize:
        records = iter_sample(records, maxsize)
    SeqIO.write(records, outfile, 'fasta')


@must_work(illegal_chars='-', min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def cluster_consensus(infiles, outfile, directory, prefix):
    """Alignment-free cluster consensus."""
    seq_id = next(SeqIO.parse(infiles[0], 'fastq')).id
    label = re.split("_[0-9]+$", seq_id)[0]
    log_ins = np.log10(globals_.config.getfloat('Parameters', 'consensus_p_ins'))
    log_del = np.log10(globals_.config.getfloat('Parameters', 'consensus_p_del'))
    kwargs = {
        'julia': globals_.config.get('Paths', 'julia'),
        'script': globals_.config.get('Paths', 'consensus_script'),
        'pattern': os.path.join(directory, "*.fastq"),
        'prefix': '{}_'.format(label),
        'outfile': outfile,
        'batch': globals_.config.get('Parameters', 'consensus_batch_size'),
        'log_ins': log_ins,
        'log_del': log_del,
        }
    cmd = ('{julia} {script} --prefix \'{prefix}\' --batch \'{batch}\''
           ' \'{log_ins}\' \'{log_del}\' \'{pattern}\' > {outfile}').format(**kwargs)
    return maybe_qsub(cmd, infiles, outfile, name="cluster-consensus-{}".format(prefix))


# making wrappers like this is necessary because nested function
# definitions are not picklable.

@must_work(seq_ids=True)
@report_wrapper
def cat_wrapper_ids(infiles, outfile):
    cat_files(infiles, outfile)


@report_wrapper
def cat_wrapper(infiles, outfile):
    cat_files(infiles, outfile)


def is_in_frame(seq):
    if len(seq) % 3 != 0:
        return False
    t = seq.translate()
    if len(t) != len(seq) / 3:
        return False
    if '*' in t:
        return False
    return True


def new_record_id(record, new_id):
    result = record[:]
    result.id = new_id
    result.name = ""
    result.description = ""
    return result


@must_work()
@report_wrapper
def inframe_hqcs(infile, outfile):
    """Filter to in-frame hqcs sequences with no stop codons"""
    records = SeqIO.parse(infile, 'fasta')
    result = list(new_record_id(r, '{}_inframe'.format(r.id))
                  for r in records if is_in_frame(r.seq))
    SeqIO.write(result, outfile, 'fasta')


def pause_for_editing_inframe_hqcs():
    infile = os.path.join(pipeline_dir, 'inframe_hqcs.fasta')
    outfile = os.path.join(pipeline_dir, 'inframe_hqcs.edited.fasta')
    if globals_.config.getboolean('Tasks', 'use_inframe_hqcs'):
        if globals_.config.getboolean('Tasks', 'pause_for_inframe_hqcs'):
            input('Paused for manual editing of {}'
                  '\nPress Enter to continue.'.format(outfile))
        # check that result is not empty
        n = sum(1 for r in SeqIO.parse(outfile, 'fasta'))
        if n == 0:
            raise Exception('{} is empty'.format(outfile))
    else:
        # just use global reference file
        dbfile = globals_.config.get('Parameters', 'reference_db')
        shutil.copyfile(dbfile, outfile)


@must_work()
@report_wrapper
def hqcs_db_search(infiles, outfile):
    infile, dbfile = infiles
    check_basename(dbfile, 'inframe_hqcs.edited.fasta')
    identity = globals_.config.get('Parameters', 'reference_identity')
    max_accepts = globals_.config.get('Parameters', 'max_accepts')
    max_rejects = globals_.config.get('Parameters', 'max_rejects')
    cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
           ' -fastapairs {outfile} -alnout {outfile}.human -userout {outfile}.calns'
           ' -userfields caln -top_hit_only -strand both'
           ' -maxaccepts {max_accepts} -maxrejects {max_rejects}')
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, db=dbfile, id=identity, outfile=outfile,
                     max_accepts=max_accepts, max_rejects=max_rejects)
    return maybe_qsub(cmd, infile, outfiles=outfile, name='hqcs-db-search')


def shift_correction_helper(infile, outfile, keep, correct, name=None):
    python = globals_.config.get('Paths', 'python')
    script = os.path.join(globals_.script_dir, 'correct_shifts.py')
    keep_option = "--keep" if keep else ""
    correct_option = "--correct-dels" if correct else ""
    calns_file = "{}.calns".format(infile)
    discard_file = "{}.discarded".format(outfile)
    summary_file = "{}.summary".format(outfile)
    cmd = ("{python} {script} {keep_option} {correct_option} --calns={calns_file}"
           " --discard={discard_file} --summary={summary_file}"
           " {infile} {outfile}")
    cmd = cmd.format(python=python, script=script, keep_option=keep_option,
                     correct_option=correct_option,
                     calns_file=calns_file, discard_file=discard_file,
                     summary_file=summary_file, infile=infile, outfile=outfile)
    return maybe_qsub(cmd, infile, outfile, name=name)

@must_work(in_frame=True, illegal_chars='-')
@report_wrapper
def hqcs_shift_correction(infile, outfile):
    return shift_correction_helper(infile, outfile, keep=False, correct=True,
                                   name="hqcs-shift-correction")


@must_work(min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def unique_hqcs(infile, outfile):
    cmd = ('{usearch} -derep_fulllength {infile} -fastaout {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    return maybe_qsub(cmd, infile, outfiles=outfile, name='unique_hqcs')


# FIXME: this is run with remote job limiter, but part of its task is run locally
@must_work()
@report_wrapper
def compute_copynumbers(infiles, outfile):
    ccsfile, hqcsfile = infiles
    # make sure this suffix changes depending on what task comes before
    check_suffix(hqcsfile, '.uniques.fasta')
    check_suffix(ccsfile, '.lenfilter.fastq')
    pairfile = '{}.pairs'.format(remove_suffix(outfile, '.tsv'))
    result = usearch_hqcs_ids(ccsfile, pairfile, hqcsfile, name='compute-copynumber')
    with open(pairfile) as f:
        pairs = list(line.strip().split("\t") for line in f.readlines())
    hqcs_counts = defaultdict(lambda: 0)
    for pair in pairs:
        if len(pair) != 2:
            warnings.warn('CCS {} did not match any HQCS'.format(pair))
            continue
        ccs_id, hqcs_id = pair
        hqcs_counts[hqcs_id] += 1
    # deal with hqcs sequences with no copynumber by giving them 0
    ids = list(r.id for r in SeqIO.parse(hqcsfile, 'fasta'))
    for i in ids:
        if i not in hqcs_counts:
            hqcs_counts[i] = 0
    with open(outfile, 'w') as handle:
        for id_, count in hqcs_counts.items():
            handle.write('{}\t{}\n'.format(id_, count))
    return result


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


def make_consensus_pipeline(name=None):
    if name is None:
        name = "consensus_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    cluster_task = pipeline.transform(cluster,
                                      input=None,
                                      filter=formatter(r'.*/(?P<LABEL>.+).qfilter'),
                                      output=os.path.join(pipeline_dir, '{LABEL[0]}.clustered.uc'))
    cluster_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    fastq_clusters_task = pipeline.subdivide(fastq_clusters,
                                             input=None,
                                             filter=formatter(r'.*/(?P<LABEL>.+).qfilter'),
                                             add_inputs=add_inputs(os.path.join(pipeline_dir, '{LABEL[0]}.clustered.uc')),
                                             output=os.path.join(pipeline_dir, '{LABEL[0]}.clusters/*.fastq'),
                                             extras=[os.path.join(pipeline_dir, '{LABEL[0]}.clusters')])
    fastq_clusters_task.jobs_limit(n_remote_jobs, remote_job_limiter)
    fastq_clusters_task.mkdir(cluster_task,
                             formatter(r'.*/(?P<LABEL>.+).clustered.uc'),
                             '{path[0]}/{LABEL[0]}.clusters')

    sample_clusters_task = pipeline.transform(sample_clusters,
                                              input=fastq_clusters_task,
                                              filter=suffix('.fastq'),
                                              output='.sampled.fastq')
    sample_clusters_task.jobs_limit(n_local_jobs, local_job_limiter)

    cluster_consensus_task = pipeline.collate(cluster_consensus,
                                              input=sample_clusters_task,
                                              filter=formatter(r'.*/(?P<LABEL>.+).clusters'),
                                              output=os.path.join(pipeline_dir, '{LABEL[0]}.consensus.fasta'),
                                              extras=['{path[0]}',
                                                      '{LABEL[0]}'])
    cluster_consensus_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    inframe_hqcs_task = pipeline.transform(inframe_hqcs,
                                           input=cluster_consensus_task,
                                           filter=suffix('.fasta'),
                                           output='.inframe.fasta')
    inframe_hqcs_task.jobs_limit(n_local_jobs, local_job_limiter)

    cat_inframe_hqcs_task = pipeline.merge(cat_wrapper_ids,
                                           input=inframe_hqcs_task,
                                           output=os.path.join(pipeline_dir, 'inframe_hqcs.fasta'))
    cat_inframe_hqcs_task.jobs_limit(n_local_jobs, local_job_limiter)

    copy_inframe_hqcs_task = pipeline.transform(copy_file,
                                                input=cat_inframe_hqcs_task,
                                                filter=suffix('.fasta'),
                                                output='.edited.fasta')
    copy_inframe_hqcs_task.jobs_limit(n_local_jobs, local_job_limiter)
    copy_inframe_hqcs_task.posttask(pause_for_editing_inframe_hqcs)

    hqcs_db_search_task = pipeline.transform(hqcs_db_search,
                                             input=cluster_consensus_task,
                                             add_inputs=add_inputs(copy_inframe_hqcs_task),
                                             filter=suffix('.fasta'),
                                             output='.refpairs.fasta')
    hqcs_db_search_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    hqcs_shift_correction_task = pipeline.transform(hqcs_shift_correction,
                                                    input=hqcs_db_search_task,
                                                    filter=suffix('.fasta'),
                                                    output='.shifted.fasta')
    hqcs_shift_correction_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    unique_hqcs_task = pipeline.transform(unique_hqcs,
                                          input=hqcs_shift_correction_task,
                                          filter=suffix('.fasta'),
                                          output='.uniques.fasta')
    unique_hqcs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    compute_copynumbers_task = pipeline.transform(compute_copynumbers,
                                                  input=None,
                                                  filter=formatter(r'.*/(?P<LABEL>.+).qfilter'),
                                                  add_inputs=add_inputs(
            os.path.join(pipeline_dir, '{LABEL[0]}.consensus.refpairs.shifted.uniques.fasta')),
                                                  output=os.path.join(pipeline_dir, '{LABEL[0]}.copynumbers.tsv'))
    compute_copynumbers_task.jobs_limit(n_remote_jobs, remote_job_limiter)
    compute_copynumbers_task.follows(unique_hqcs_task)

    merge_copynumbers_task = pipeline.merge(cat_wrapper,
                                            name='cat_copynumbers',
                                            input=compute_copynumbers_task,
                                            output=os.path.join(pipeline_dir, 'copynumbers.tsv'))
    merge_copynumbers_task.jobs_limit(n_local_jobs, local_job_limiter)

    cat_all_hqcs_task = pipeline.merge(cat_all_hqcs,
                                       input=unique_hqcs_task,
                                       output=os.path.join(pipeline_dir, "hqcs.fasta"))
    cat_all_hqcs_task.jobs_limit(n_local_jobs, local_job_limiter)

    pipeline.set_head_tasks([cluster_task, fastq_clusters_task, compute_copynumbers_task])
    pipeline.set_tail_tasks([cat_all_hqcs_task, merge_copynumbers_task])

    for task in pipeline.head_tasks:
        task.mkdir(pipeline_dir)

    return pipeline
