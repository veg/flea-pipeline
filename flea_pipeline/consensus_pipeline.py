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


pipeline_dir = os.path.join(globals_.results_dir, "consensus")


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


@must_work(maybe=True)
@report_wrapper
def mafft(infile, outfile):
    stderr = '{}.stderr'.format(outfile)
    cmd = '/usr/local/bin/mafft-fftns --ep 0.5 --quiet --preservecase {}'.format(infile)
    result = run_command(cmd, infile, outfile, stdout=outfile,
                         stderr=stderr, name="mafft")
    return result


@must_work(maybe=True, illegal_chars='-')
@report_wrapper
def alignment_consensus(infile, outfile):
    n = re.search("cluster_([0-9]+)", infile).group(1)
    seq_id = next(SeqIO.parse(infile, 'fasta')).id
    label = re.split("_[0-9]+$", seq_id)[0]
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "DNAcons.py"),
        'infile': infile,
        'outfile': outfile,
        'ambifile': '{}.info'.format(outfile[:-len('.fasta')]),
        'id_str': "{label}_hqcs_{n}".format(label=label, n=n),
        }
    cmd = ("{python} {script} -o {outfile} --ambifile {ambifile}"
           " --id {id_str} {infile}".format(**kwargs))
    return run_command(cmd, infile, outfile, name="cluster-consensus")


def shift_correction_helper(infile, outfile, keep, correct, name=None):
    if name is None:
        name = 'shift_correction_helper'
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
    return run_command(cmd, infile, outfile, name=name)


@must_work(in_frame=True, illegal_chars='-')
@report_wrapper
def shift_correction(infile, outfile):
    return shift_correction_helper(infile, outfile, keep=False, correct=True,
                                   name="shift-correction")


def fastq_clusters_helper(infiles, outfiles, outdir, pattern, fasta=False):
    for f in outfiles:
        os.unlink(f)
    ucfile, fastqfile = infiles
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "cluster_fastq.py"),
        'ucfile': ucfile,
        'fastqfile': fastqfile,
        'outdir': outdir,
        'fasta': '--fasta' if fasta else '',
        'minsize': globals_.config.get('Parameters', 'min_cluster_size'),
        }
    cmd = ("{python} {script} {fasta} --minsize {minsize}"
           " {ucfile} {fastqfile} {outdir}".format(**kwargs))
    return run_command(cmd, infiles, [], name="cluster-fastq")


@must_work(many=True)
@report_wrapper
def fastq_clusters(infiles, outfiles, outdir, pattern):
    fastq_clusters_helper(infiles, outfiles, outdir, pattern)


@must_work(many=True)
@report_wrapper
def fasta_clusters(infiles, outfiles, outdir, pattern):
    fastq_clusters_helper(infiles, outfiles, outdir, pattern, fasta=True)


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

def get_format(f):
    return os.path.splitext(f)[1][1:]

@must_work()
@report_wrapper
def sample_clusters(infile, outfile):
    # TODO: sample by quality score
    format = get_format(infile)
    n = sum(1 for r in SeqIO.parse(infile, format))
    records = SeqIO.parse(infile, format)
    maxsize = int(globals_.config.get('Parameters', 'max_cluster_size'))
    if n > maxsize:
        records = iter_sample(records, maxsize)
    SeqIO.write(records, outfile, format)


@must_work(illegal_chars='-', min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def cluster_consensus(infiles, outfile, directory, prefix):
    """Alignment-free cluster consensus."""
    ppn = 1 if globals_.run_locally else globals_.ppn
    options = ''
    if ppn > 1 and globals_.config.getboolean('Parameters', 'consensus_multiprocess'):
        options = '-p {}'.format(ppn)
    indel_file = '{}.indel-probs.txt'.format(remove_suffix(outfile, '.fastq'))
    kwargs = {
        'julia': globals_.config.get('Paths', 'julia'),
        'script': globals_.config.get('Paths', 'consensus_script'),
        'options': options,
        'phred_cap': globals_.config.get('Parameters', 'phred_cap'),
        'maxiters': globals_.config.get('Parameters', 'consensus_max_iters'),
        'seq_errors': globals_.config.get('Parameters', 'seq_errors'),
        'indel_file': indel_file,
        'pattern': os.path.join(directory, "*.sampled.fastq"),
        'outfile': outfile,
        }
    cmd = ('{julia} {options} {script} '
           ' --phred-cap \'{phred_cap}\''
           ' --max-iters \'{maxiters}\''
           ' \'{seq_errors}\' \'{pattern}\' {outfile}').format(**kwargs)
    return run_command(cmd, infiles, outfile,
                       ppn=ppn, name="cluster-consensus-{}".format(prefix))


@must_work(illegal_chars='-', min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def cluster_consensus_with_ref(infiles, outfile, directory, prefix):
    """Alignment-free cluster consensus."""
    # TODO: if reference is this cluster's HQCS, done
    reffile = output=os.path.join(pipeline_dir, "hqcs_refs_inframe_db.fasta")
    refmapfile = os.path.join(pipeline_dir,
                              "{}.hqcs-ref-id-map.txt".format(prefix))
    seq_id = next(SeqIO.parse(infiles[0], 'fastq')).id
    label = re.split("_[0-9]+$", seq_id)[0]
    ppn = 1 if globals_.run_locally else globals_.ppn
    options = ''
    if ppn > 1 and globals_.config.getboolean('Parameters', 'consensus_multiprocess'):
        options = '-p {}'.format(ppn)
    kwargs = {
        'julia': globals_.config.get('Paths', 'julia'),
        'script': globals_.config.get('Paths', 'consensus_script'),
        'options': options,
        'prefix': '{}_hqcs_'.format(label),
        'phred_cap': globals_.config.get('Parameters', 'phred_cap'),
        'maxiters': globals_.config.get('Parameters', 'consensus_max_iters'),
        'reffile': reffile,
        'refmapfile': refmapfile,
        'ref_errors': globals_.config.get('Parameters', 'ref_errors'),
        'seq_errors': globals_.config.get('Parameters', 'seq_errors'),
        'pattern': os.path.join(directory, "*.sampled.fastq"),
        'outfile': outfile,
        }
    cmd = ('{julia} {options} {script} '
           ' --prefix \'{prefix}\''
           ' --keep-unique-name'
           ' --phred-cap \'{phred_cap}\''
           ' --max-iters \'{maxiters}\''
           ' --reference \'{reffile}\''
           ' --reference-map \'{refmapfile}\''
           ' --ref-errors \'{ref_errors}\''
           ' \'{seq_errors}\' \'{pattern}\' {outfile}').format(**kwargs)
    return run_command(cmd, infiles, outfile, ppn=ppn, name="cluster-consensus-{}".format(prefix))


@must_work()
@report_wrapper
def all_hq_hqcs(infile, outfile):
    max_error_rate = globals_.config.getfloat('Parameters', 'hqcs_max_err_rate')
    max_base_error_rate = globals_.config.getfloat('Parameters', 'hqcs_max_base_err_rate')
    records = SeqIO.parse(infile, 'fastq')
    # FIXME: what if none are left after filtering???
    def filt(r):
        phreds = np.array(r.letter_annotations["phred_quality"])
        errors = 10 ** (-phreds / 10.0)
        return sum(errors) <= max_error_rate and max(errors) <= max_base_error_rate
    to_keep = (r for r in records if filt(r))
    SeqIO.write(to_keep, outfile, 'fastq')

# making wrappers like this is necessary because nested function
# definitions are not picklable.

@must_work(seq_ids=True)
@report_wrapper
def cat_wrapper_ids(infiles, outfile):
    cat_files(infiles, outfile)


@report_wrapper
def cat_wrapper(infiles, outfile):
    cat_files(infiles, outfile)


def is_in_frame(seq, allow_stop_codons):
    if len(seq) % 3 != 0:
        return False
    t = seq.translate()
    if len(t) != len(seq) / 3:
        return False
    if '*' in t and not allow_stop_codons:
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
def inframe_nostops(infile, outfile):
    """Filter to in-frame hqcs sequences with no stop codons"""
    records = SeqIO.parse(infile, 'fastq')
    result = list(new_record_id(r, '{}_inframe'.format(r.id))
                  for r in records if is_in_frame(r.seq, False))
    SeqIO.write(result, outfile, 'fasta')


@must_work()
@report_wrapper
def inframe_allow_stops(infile, outfile):
    """Filter to in-frame hqcs sequences, allowing stop codons"""
    format = get_format(infile)
    records = SeqIO.parse(infile, format)
    result = list(r for r in records if is_in_frame(r.seq, True))
    SeqIO.write(result, outfile, 'fasta')


@report_wrapper
def outframe_hqcs(infile, outfile):
    """Filter to sequences that are not in frame"""
    format = get_format(infile)
    records = SeqIO.parse(infile, format)
    result = list(r for r in records if not is_in_frame(r.seq, False))
    SeqIO.write(result, outfile, 'fasta')


#@must_work(min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@must_work()
@report_wrapper
def unique_hqcs(infile, outfile):
    cmd = ('{usearch} -derep_fulllength {infile} -fastaout {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    return run_command(cmd, infile, outfiles=outfile, name='unique_hqcs')


def pause_for_editing_inframe_hqcs():
    infile = os.path.join(pipeline_dir, 'hqcs_refs_inframe_db.fasta')
    outfile = os.path.join(pipeline_dir, 'hqcs_refs_inframe_db.edited.fasta')
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

def hqcs_db_search_helper(infiles, outfile, fasta=False):
    infile, dbfile = infiles
    check_basename(dbfile, 'hqcs_refs_inframe_db(.edited)?.fasta')
    identity = globals_.config.get('Parameters', 'reference_identity')
    max_accepts = globals_.config.get('Parameters', 'max_accepts')
    max_rejects = globals_.config.get('Parameters', 'max_rejects')
    if fasta:
        output_args = ('-fastapairs {outfile} -userout {outfile}.calns '
                       ' -userfields caln'.format(outfile=outfile))
    else:
        output_args = '-userout {} -userfields query+target'.format(outfile)
    cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
           ' {output_args} -top_hit_only -strand both'
           ' -maxaccepts {max_accepts} -maxrejects {max_rejects}')
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, db=dbfile, id=identity,
                     output_args=output_args,
                     max_accepts=max_accepts, max_rejects=max_rejects)
    return run_command(cmd, infile, outfiles=outfile, name='hqcs-db-search')


@must_work()
@report_wrapper
def hqcs_db_search_fasta(infiles, outfile):
    return hqcs_db_search_helper(infiles, outfile, fasta=True)


@must_work()
@report_wrapper
def hqcs_db_search_ids(infiles, outfile):
    return hqcs_db_search_helper(infiles, outfile, fasta=False)


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

    use_rifraf = globals_.config.getboolean('Parameters', 'use_rifraf')

    if use_rifraf:
        # TODO: use subpipelines for all of these
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

        # next few tasks: find cluster consensus without reference
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

        all_hq_hqcs_task = pipeline.transform(all_hq_hqcs,
                                               input=cluster_consensus_task,
                                               filter=suffix('.fastq'),
                                               output='.allhq.fastq')

        inframe_refhqcs_task = pipeline.transform(inframe_nostops,
                                                  input=all_hq_hqcs_task,
                                                  filter=suffix('.fastq'),
                                                  output='.inframe.fasta')

        unique_hqcs_ref_task = pipeline.transform(unique_hqcs,
                                                  name="unique_hqcs_ref",
                                                  input=inframe_refhqcs_task,
                                                  filter=suffix('.fasta'),
                                                  output='.uniques.fasta')

        cat_all_hqcs_ref_task = pipeline.merge(cat_all_hqcs,
                                               name="cat_all_hqcs_ref",
                                               input=unique_hqcs_ref_task,
                                               output=os.path.join(pipeline_dir, "hqcs_refs_inframe_db.fasta"))

        # next few tasks: find cluster consensus with reference, using
        # previously-found consensus as reference
        hqcs_db_search_task = pipeline.transform(hqcs_db_search_ids,
                                                 input=cluster_consensus_task,
                                                 filter=suffix('.refhqcs.fastq'),
                                                 output=".hqcs-ref-id-map.txt",
                                                 add_inputs=add_inputs(cat_all_hqcs_ref_task))

        cluster_consensus_with_ref_task = pipeline.collate(cluster_consensus_with_ref,
                                                           input=sample_clusters_task,
                                                           filter=formatter(r'(.*)/(?P<NAME>[a-zA-Z0-9-_]*)\.clusters(.*)'),
                                                           output=os.path.join(pipeline_dir, '{NAME[0]}.hqcs.fastq'),
                                                           extras=['{path[0]}',
                                                                   '{NAME[0]}'])
        cluster_consensus_with_ref_task.follows(hqcs_db_search_task)

        previous_task = cluster_consensus_with_ref
    else:
        # use MAFFT and shift corection
        fasta_clusters_task = pipeline.subdivide(fasta_clusters,
                                                 input=cluster_task,
                                                 filter=formatter(basename_regex),
                                                 add_inputs=add_inputs(ccs_pattern),
                                                 output=os.path.join(pipeline_dir, '{NAME[0]}.clusters/*.raw.fasta'),
                                                 extras=[os.path.join(pipeline_dir, '{NAME[0]}.clusters'),
                                                         'cluster_[0-9]+.raw.fasta'])
        fasta_clusters_task.mkdir(cluster_task,
                                 formatter(r'.*/(?P<NAME>.+).clustered.uc'),
                                 '{path[0]}/{NAME[0]}.clusters')

        sample_clusters_task = pipeline.transform(sample_clusters,
                                                  input=fasta_clusters_task,
                                                  filter=suffix('.raw.fasta'),
                                                  output='.sampled.fasta')

        align_clusters_task = pipeline.transform(mafft,
                                                 name="align_clusters",
                                                 input=sample_clusters_task,
                                                 filter=suffix('.fasta'),
                                                 output='.aligned.fasta')

        alignment_consensus_task = pipeline.transform(alignment_consensus,
                                                    input=align_clusters_task,
                                                    filter=suffix('.fasta'),
                                                    output='.hqcs.fasta')

        cat_clusters_task = pipeline.collate(cat_wrapper_ids,
                                             name="cat_clusters",
                                             input=alignment_consensus_task,
                                             filter=formatter(),
                                             output='{path[0]}.hqcs.fasta')

        if globals_.config.getboolean('Parameters', 'do_shift_correction'):
            inframe_hqcs_refs_task = pipeline.transform(inframe_nostops,
                                                        name="inframe_hqcs_refs",
                                                       input=cat_clusters_task,
                                                       filter=suffix('.fasta'),
                                                       output='.inframe.fasta')

            cat_inframe_hqcs_refs_task = pipeline.merge(cat_wrapper_ids,
                                                        input=inframe_hqcs_refs_task,
                                                        output=os.path.join(pipeline_dir,
                                                                            'hqcs_refs_inframe_db.fasta'))

            copy_inframe_hqcs_task = pipeline.transform(copy_file,
                                                        input=cat_inframe_hqcs_refs_task,
                                                        filter=suffix('.fasta'),
                                                        output='.edited.fasta')
            copy_inframe_hqcs_task.posttask(pause_for_editing_inframe_hqcs)

            hqcs_db_search_task = pipeline.transform(hqcs_db_search_fasta,
                                                     input=cat_clusters_task,
                                                     add_inputs=add_inputs(copy_inframe_hqcs_task),
                                                     filter=suffix('.fasta'),
                                                     output='.refpairs.fasta')

            shift_correction_task = pipeline.transform(shift_correction,
                                                       input=hqcs_db_search_task,
                                                       filter=suffix('.fasta'),
                                                       output='.shifted.fasta')

            previous_task = shift_correction_task
        else:
            previous_task = cat_clusters_task
    # end

    inframe_hqcs_task = pipeline.transform(inframe_allow_stops,
                                           input=previous_task,
                                           filter=suffix('.fastq' if use_rifraf else '.fasta'),
                                           output='.inframe.fasta')

    outframe_hqcs_task = pipeline.transform(outframe_hqcs,
                                            input=previous_task,
                                            filter=suffix('.fastq' if use_rifraf else '.fasta'),
                                            output='.outframe.fasta')

    unique_hqcs_task = pipeline.transform(unique_hqcs,
                                          input=inframe_hqcs_task,
                                          filter=suffix('.fasta'),
                                          output='.uniques.fasta')

    cat_all_hqcs_task = pipeline.merge(cat_all_hqcs,
                                       name="cat_all_hqcs",
                                       input=unique_hqcs_task,
                                       output=os.path.join(pipeline_dir, "hqcs.fasta"))

    compute_copynumbers_task = pipeline.transform(compute_copynumbers,
                                                  input=unique_hqcs_task,
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
