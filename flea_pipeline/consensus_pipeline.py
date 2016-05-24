import os
import re
import random
from collections import defaultdict
import warnings
import shutil

import numpy as np
from ruffus import Pipeline, suffix, formatter, add_inputs
from Bio import SeqIO

from flea_pipeline.util import run_command, cat_files
from flea_pipeline.util import must_work, report_wrapper
from flea_pipeline.util import check_suffix, check_basename
from flea_pipeline.util import remove_suffix
from flea_pipeline.util import usearch_hqcs_ids

import flea_pipeline.pipeline_globals as globals_


pipeline_dir = os.path.join(globals_.data_dir, "consensus")


@must_work()
@report_wrapper
def make_input(infile, outfile):
    print(infile, outfile)
    if os.path.exists(outfile):
        os.unlink(outfile)
    os.symlink(infile, outfile)


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
    return run_command(cmd, infile, outfile, name="cluster")


@must_work(many=True)
@report_wrapper
def fastq_clusters(infiles, outfiles, outdir, pattern):
    for f in outfiles:
        os.unlink(f)
    ucfile, fastqfile = infiles
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
    return run_command(cmd, infiles, [], name="cluster-fastq")


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
    SeqIO.write(records, outfile, 'fastq')


@must_work(illegal_chars='-', min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def cluster_consensus(infiles, outfile, directory, prefix):
    """Alignment-free cluster consensus."""
    seq_id = next(SeqIO.parse(infiles[0], 'fastq')).id
    label = re.split("_[0-9]+$", seq_id)[0]
    log_ins = np.log10(globals_.config.getfloat('Parameters', 'consensus_p_ins'))
    log_del = np.log10(globals_.config.getfloat('Parameters', 'consensus_p_del'))
    ppn = 1 if globals_.run_locally else globals_.ppn
    options = ''
    if ppn > 1:
        options = '-p {}'.format(ppn)
    kwargs = {
        'julia': globals_.config.get('Paths', 'julia'),
        'script': globals_.config.get('Paths', 'consensus_script'),
        'options': options,
        'pattern': os.path.join(directory, "*.sampled.fastq"),
        'refpattern': globals_.config.get('Parameters', 'reference_dna'),
        'prefix': '{}_hqcs_'.format(label),
        'outfile': outfile,
        'batch': globals_.config.get('Parameters', 'consensus_batch_size'),
        'log_ins': log_ins,
        'log_del': log_del,
        }
    cmd = ('{julia} {options} {script} --prefix \'{prefix}\''
           ' --reference \'{refpattern}\''
           ' --keep-unique-name --batch \'{batch}\''
           ' \'{log_ins}\' \'{log_del}\' \'{pattern}\' > {outfile}').format(**kwargs)
    return run_command(cmd, infiles, outfile, ppn=ppn, name="cluster-consensus-{}".format(prefix))


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
    return run_command(cmd, infile, outfiles=outfile, name='hqcs-db-search')


@must_work(min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def unique_hqcs(infile, outfile):
    cmd = ('{usearch} -derep_fulllength {infile} -fastaout {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    return run_command(cmd, infile, outfiles=outfile, name='unique_hqcs')


# FIXME: this is run with remote job limiter, but part of its task is run locally
@must_work()
@report_wrapper
def compute_copynumbers(infiles, outfile):
    hqcsfile, ccsfile = infiles
    # make sure this suffix changes depending on what task comes before
    check_suffix(hqcsfile, '.uniques.fasta')
    check_suffix(ccsfile, '.ccs.fastq')
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

    inputs_regex = r'.*/.*(?P<LABEL>{}).*'.format('|'.join(t.label for t in globals_.timepoints))
    regex = r'.*/(?P<LABEL>{}).*'.format('|'.join(t.label for t in globals_.timepoints))
    ccs_pattern = os.path.join(pipeline_dir, '{LABEL[0]}.ccs.fastq')

    make_inputs_task = pipeline.transform(make_input,
                                          input=None,
                                          filter=formatter(inputs_regex),
                                          output=os.path.join(pipeline_dir, '{LABEL[0]}.ccs.fastq'))
    make_inputs_task.mkdir(pipeline_dir)

    cluster_task = pipeline.transform(cluster,
                                      input=make_inputs_task,
                                      filter=formatter(regex),
                                      output=os.path.join(pipeline_dir, '{LABEL[0]}.clustered.uc'))

    fastq_clusters_task = pipeline.subdivide(fastq_clusters,
                                             input=cluster_task,
                                             filter=formatter(regex),
                                             add_inputs=add_inputs(ccs_pattern),
                                             output=os.path.join(pipeline_dir, '{LABEL[0]}.clusters/*.raw.fastq'),
                                             extras=[os.path.join(pipeline_dir, '{LABEL[0]}.clusters'),
                                                     'cluster_[0-9]+.raw.fastq'])
    fastq_clusters_task.mkdir(cluster_task,
                             formatter(r'.*/(?P<LABEL>.+).clustered.uc'),
                             '{path[0]}/{LABEL[0]}.clusters')

    sample_clusters_task = pipeline.transform(sample_clusters,
                                              input=fastq_clusters_task,
                                              filter=suffix('.raw.fastq'),
                                              output='.sampled.fastq')

    cluster_consensus_task = pipeline.collate(cluster_consensus,
                                              input=sample_clusters_task,
                                              filter=formatter(regex),
                                              output=os.path.join(pipeline_dir, '{LABEL[0]}.hqcs.fasta'),
                                              extras=['{path[0]}',
                                                      '{LABEL[0]}'])

    inframe_hqcs_task = pipeline.transform(inframe_hqcs,
                                           input=cluster_consensus_task,
                                           filter=suffix('.fasta'),
                                           output='.inframe.fasta')

    # TODO: redo consensus with orthologous references

    unique_hqcs_task = pipeline.transform(unique_hqcs,
                                          input=inframe_hqcs_task,
                                          filter=suffix('.fasta'),
                                          output='.uniques.fasta')

    compute_copynumbers_task = pipeline.transform(compute_copynumbers,
                                                  input=unique_hqcs_task,
                                                  filter=formatter(regex),
                                                  add_inputs=add_inputs(ccs_pattern),
                                                  output=os.path.join(pipeline_dir, '{LABEL[0]}.copynumbers.tsv'))

    merge_copynumbers_task = pipeline.merge(cat_wrapper,
                                            name='cat_copynumbers',
                                            input=compute_copynumbers_task,
                                            output=os.path.join(pipeline_dir, 'copynumbers.tsv'))

    cat_all_hqcs_task = pipeline.merge(cat_all_hqcs,
                                       input=unique_hqcs_task,
                                       output=os.path.join(pipeline_dir, "hqcs.fasta"))

    pipeline.set_head_tasks([make_inputs_task])
    pipeline.set_tail_tasks([cat_all_hqcs_task, merge_copynumbers_task])

    return pipeline
