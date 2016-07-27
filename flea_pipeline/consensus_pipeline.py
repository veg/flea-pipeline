import os
import re
import random
from collections import defaultdict
import warnings
import shutil
import csv

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
    ppn = 1 if globals_.run_locally else globals_.ppn
    options = ''
    if ppn > 1:
        options = '-p {}'.format(ppn)
    indel_file = '{}.indel-probs.txt'.format(remove_suffix(outfile, '.fastq'))
    kwargs = {
        'julia': globals_.config.get('Paths', 'julia'),
        'script': globals_.config.get('Paths', 'consensus_script'),
        'options': options,
        'batch': globals_.config.get('Parameters', 'consensus_batch_size'),
        'maxiters': globals_.config.get('Parameters', 'consensus_max_iters'),
        'indel_file': indel_file,
        'pattern': os.path.join(directory, "*.sampled.fastq"),
        'outfile': outfile,
        }
    cmd = ('{julia} {options} {script} --batch \'{batch}\''
           ' --max-iters \'{maxiters}\''
           ' --indel-file \'{indel_file}\''
           ' \'{pattern}\' {outfile}').format(**kwargs)
    return run_command(cmd, infiles, outfile,
                       ppn=ppn, name="cluster-consensus-{}".format(prefix))


@must_work(illegal_chars='-', min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def cluster_consensus_with_ref(infiles, outfile, directory, prefix):
    """Alignment-free cluster consensus."""
    # TODO: if reference is this cluster's HQCS, done
    reffile = output=os.path.join(pipeline_dir, "inframe-db.fastq")
    refmapfile = os.path.join(pipeline_dir,
                              "{}.hqcs-ref-id.txt".format(prefix))
    seq_id = next(SeqIO.parse(infiles[0], 'fastq')).id
    label = re.split("_[0-9]+$", seq_id)[0]
    ppn = 1 if globals_.run_locally else globals_.ppn
    options = ''
    if ppn > 1:
        options = '-p {}'.format(ppn)
    indel_file = '{}.indel-probs.txt'.format(remove_suffix(outfile, '.fastq'))
    kwargs = {
        'julia': globals_.config.get('Paths', 'julia'),
        'script': globals_.config.get('Paths', 'consensus_script'),
        'options': options,
        'prefix': '{}_hqcs_'.format(label),
        'batch': globals_.config.get('Parameters', 'consensus_batch_size'),
        'maxiters': globals_.config.get('Parameters', 'consensus_max_iters'),
        'indel_file': indel_file,
        'reffile': reffile,
        'refmapfile': refmapfile,
        'pattern': os.path.join(directory, "*.sampled.fastq"),
        'outfile': outfile,
        }
    cmd = ('{julia} {options} {script} --prefix \'{prefix}\''
           ' --keep-unique-name --batch \'{batch}\''
           ' --max-iters \'{maxiters}\''
           ' --indel-file \'{indel_file}\''
           ' --reference \'{reffile}\''
           ' --reference-map \'{refmapfile}\''
           ' \'{pattern}\' {outfile}').format(**kwargs)
    return run_command(cmd, infiles, outfile, ppn=ppn, name="cluster-consensus-{}".format(prefix))


@must_work()
@report_wrapper
def no_indel_hqcs(infile, outfile):
    indel_file = '{}.indel-probs.txt'.format(remove_suffix(infile, '.fastq'))
    max_error_rate = globals_.config.getfloat('Parameters', 'hqcs_max_err_rate')
    max_base_error_rate = globals_.config.getfloat('Parameters', 'hqcs_max_base_err_rate')
    to_keep = set()
    with open(indel_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if (float(row["max_indel_p"]) <= max_base_error_rate and
                float(row["indel_rate"]) <= max_base_error_rate):
                to_keep.add(row["name"])
    records = (r for r in SeqIO.parse(infile, 'fastq')
               if r.id in to_keep)
    SeqIO.write(records, outfile, 'fastq')


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
    records = SeqIO.parse(infile, 'fastq')
    result = list(new_record_id(r, '{}_inframe'.format(r.id))
                  for r in records if is_in_frame(r.seq))
    SeqIO.write(result, outfile, 'fastq')


#@must_work(min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@must_work()
@report_wrapper
def unique_hqcs(infile, outfile):
    cmd = ('{usearch} -derep_fulllength {infile} -fastqout {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    return run_command(cmd, infile, outfiles=outfile, name='unique_hqcs')


# def pause_for_editing_inframe_hqcs():
#     infile = os.path.join(pipeline_dir, 'inframe_hqcs.fasta')
#     outfile = os.path.join(pipeline_dir, 'inframe_hqcs.edited.fasta')
#     if globals_.config.getboolean('Tasks', 'use_inframe_hqcs'):
#         if globals_.config.getboolean('Tasks', 'pause_for_inframe_hqcs'):
#             input('Paused for manual editing of {}'
#                   '\nPress Enter to continue.'.format(outfile))
#         # check that result is not empty
#         n = sum(1 for r in SeqIO.parse(outfile, 'fasta'))
#         if n == 0:
#             raise Exception('{} is empty'.format(outfile))
#     else:
#         # just use global reference file
#         dbfile = globals_.config.get('Parameters', 'reference_db')
#         shutil.copyfile(dbfile, outfile)


@must_work()
@report_wrapper
def hqcs_db_search(infiles, outfile):
    infile, dbfile = infiles
    check_basename(dbfile, 'inframe-db.fasta')
    identity = globals_.config.get('Parameters', 'reference_identity')
    max_accepts = globals_.config.get('Parameters', 'max_accepts')
    max_rejects = globals_.config.get('Parameters', 'max_rejects')
    cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
           ' -userout {outfile} -userfields query+target -top_hit_only -strand both'
           ' -maxaccepts {max_accepts} -maxrejects {max_rejects}')
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, db=dbfile, id=identity, outfile=outfile,
                     max_accepts=max_accepts, max_rejects=max_rejects)
    return run_command(cmd, infile, outfiles=outfile, name='hqcs-db-search')


# FIXME: this is run with remote job limiter, but part of its task is run locally
@must_work()
@report_wrapper
def compute_copynumbers(infiles, outfile):
    hqcsfile, ccsfile = infiles
    # make sure this suffix changes depending on what task comes before
    check_suffix(hqcsfile, '.uniques.fasta')
    check_suffix(ccsfile, '.ccs.fastq')
    pairfile = '{}.pairs'.format(remove_suffix(outfile, '.tsv'))
    identity = globals_.config.get('Parameters', 'copynumber_identity')
    maxqt = globals_.config.get('Parameters', 'cn_max_length_ratio')
    result = usearch_hqcs_ids(ccsfile, pairfile, hqcsfile,
                              identity, maxqt, name='compute-copynumber')
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


@must_work(seq_ids=True)
@report_wrapper
def fastq_to_fasta(infile, outfile):
    records = SeqIO.parse(infile, 'fastq')
    SeqIO.write(records, outfile, 'fasta')


@must_work()
@report_wrapper
def copy_file(infile, outfile):
    shutil.copyfile(infile, outfile)


def make_consensus_pipeline(name=None):
    if name is None:
        name = "consensus_pipeline"
    pipeline = Pipeline(name)

    basename_regex = r'(.*)/(?P<NAME>[a-zA-Z0-9-_]*)\.(.*)'
    ccs_pattern = os.path.join(pipeline_dir, '{NAME[0]}.ccs.fastq')

    make_inputs_task = pipeline.transform(make_input,
                                          input=None,
                                          filter=formatter(basename_regex),
                                          output=os.path.join(pipeline_dir, '{NAME[0]}.ccs.fastq'))
    make_inputs_task.mkdir(pipeline_dir)

    cluster_task = pipeline.transform(cluster,
                                      input=make_inputs_task,
                                      filter=formatter(basename_regex),
                                      output=os.path.join(pipeline_dir, '{NAME[0]}.clustered.uc'))

    fastq_clusters_task = pipeline.subdivide(fastq_clusters,
                                             input=cluster_task,
                                             filter=formatter(basename_regex),
                                             add_inputs=add_inputs(ccs_pattern),
                                             output=os.path.join(pipeline_dir, '{NAME[0]}.clusters/*.raw.fastq'),
                                             extras=[os.path.join(pipeline_dir, '{NAME[0]}.clusters'),
                                                     'cluster_[0-9]+.raw.fastq'])
    fastq_clusters_task.mkdir(cluster_task,
                             formatter(r'.*/(?P<NAME>.+).clustered.uc'),
                             '{path[0]}/{NAME[0]}.clusters')

    sample_clusters_task = pipeline.transform(sample_clusters,
                                              input=fastq_clusters_task,
                                              filter=suffix('.raw.fastq'),
                                              output='.sampled.fastq')

    cluster_consensus_task = pipeline.collate(cluster_consensus,
                                              input=sample_clusters_task,
                                              filter=formatter(r'(.*)/(?P<NAME>[a-zA-Z0-9-_]*)\.clusters(.*)'),
                                              output=os.path.join(pipeline_dir, '{NAME[0]}.refhqcs.fastq'),
                                              extras=['{path[0]}',
                                                      '{NAME[0]}'])

    no_indel_hqcs_task = pipeline.transform(no_indel_hqcs,
                                           input=cluster_consensus_task,
                                           filter=suffix('.fastq'),
                                           output='.no-indels.fastq')

    inframe_hqcs_task = pipeline.transform(inframe_hqcs,
                                           input=no_indel_hqcs_task,
                                           filter=suffix('.fastq'),
                                           output='.inframe.fastq')

    unique_hqcs_ref_task = pipeline.transform(unique_hqcs,
                                              name="unique_hqcs_ref",
                                              input=inframe_hqcs_task,
                                              filter=suffix('.fastq'),
                                              output='.uniques.fastq')

    cat_all_hqcs_ref_task = pipeline.merge(cat_all_hqcs,
                                           name="cat_all_hqcs_ref",
                                           input=unique_hqcs_ref_task,
                                           output=os.path.join(pipeline_dir, "inframe-db.fastq"))

    hqcs_ref_fasta_task = pipeline.transform(fastq_to_fasta,
                                             name="hqcs_ref_fasta_task",
                                             input=cat_all_hqcs_ref_task,
                                             filter=suffix('.fastq'),
                                             output='.fasta')

    hqcs_db_search_task = pipeline.transform(hqcs_db_search,
                                             input=cluster_consensus_task,
                                             filter=suffix('.refhqcs.fastq'),
                                             output=".hqcs-ref-id.txt",
                                             add_inputs=add_inputs(hqcs_ref_fasta_task))

    # 2nd round of consensus algorithm
    cluster_consensus_ref_task = pipeline.collate(cluster_consensus_with_ref,
                                                  input=sample_clusters_task,
                                                  filter=formatter(r'(.*)/(?P<NAME>[a-zA-Z0-9-_]*)\.clusters(.*)'),
                                                  output=os.path.join(pipeline_dir, '{NAME[0]}.hqcs.fastq'),
                                                  extras=['{path[0]}',
                                                          '{NAME[0]}'])
    cluster_consensus_ref_task.follows(hqcs_db_search_task)

    unique_hqcs_task = pipeline.transform(unique_hqcs,
                                          input=cluster_consensus_ref_task,
                                          filter=suffix('.fastq'),
                                          output='.uniques.fastq')

    hqcs_fasta_task = pipeline.transform(fastq_to_fasta,
                                         name="hqcs_fasta_task",
                                         input=unique_hqcs_task,
                                         filter=suffix('.fastq'),
                                         output='.fasta')

    cat_all_hqcs_task = pipeline.merge(cat_all_hqcs,
                                       name="cat_all_hqcs",
                                       input=hqcs_fasta_task,
                                       output=os.path.join(pipeline_dir, "hqcs.fasta"))

    compute_copynumbers_task = pipeline.transform(compute_copynumbers,
                                                  input=hqcs_fasta_task,
                                                  filter=formatter(basename_regex),
                                                  add_inputs=add_inputs(ccs_pattern),
                                                  output=os.path.join(pipeline_dir, '{NAME[0]}.copynumbers.tsv'))

    merge_copynumbers_task = pipeline.merge(cat_wrapper,
                                            name='cat_copynumbers',
                                            input=compute_copynumbers_task,
                                            output=os.path.join(pipeline_dir, 'copynumbers.tsv'))

    pipeline.set_head_tasks([make_inputs_task])
    pipeline.set_tail_tasks([cat_all_hqcs_task, merge_copynumbers_task])

    return pipeline
