"""HyPhy-free analysis. Currently only generates trees and sequences."""

import os
import json
from itertools import islice
from datetime import datetime
from collections import defaultdict

from ruffus import Pipeline, formatter, suffix

from Bio import SeqIO

from translate import translate
import pipeline_globals as globals_
from util import must_work, maybe_qsub, call, n_jobs, local_job_limiter, remote_job_limiter, read_single_record, cat_files
from DNAcons import consfile
from alignment_pipeline import mafft


pipeline_dir = os.path.join(globals_.data_dir, "fasttree")


@must_work()
def gapped_translate_wrapper(infiles, outfile):
    # hack to use either single infile or first of list
    if isinstance(infiles, str):
        infile = infiles
    else:
        infile, _ = infiles
    translate(infile, outfile, gapped=True)


@must_work()
def run_fasttree(infiles, outfile):
    infile, _ = infiles
    binary = globals_.config.get('Paths', 'FastTree')
    stderr = '{}.stderr'.format(outfile)
    cmd = '{} -gtr -nt {} > {} 2>{}'.format(binary, infile, outfile, stderr)
    maybe_qsub(cmd, outfiles=outfile, stdout='/dev/null', stderr='/dev/null')


def mrca(infile, recordfile, copynumber_file, outfile, oldest_id):
    """Writes records from `infile` to `recordfile` that have an ID
    corresponding to `oldest_id`. Then runs DNAcons, writing result to
    `outfile`.

    """
    len_id = len(oldest_id)
    records = SeqIO.parse(infile, "fasta")
    oldest_records = (r for r in records if r.id.startswith(oldest_id))
    SeqIO.write(oldest_records, recordfile, "fasta")
    consfile(recordfile, outfile, copynumber_file=copynumber_file,
             id_str="mrca", ungap=False,)


@must_work(illegal_chars='')
def compute_mrca(infiles, outfile):
    alignment_file, copynumber_file = infiles
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(globals_.timepoints, key=strptime)
    oldest_records_filename = '.'.join([alignment_file, "oldest_{}".format(oldest_timepoint.id)])
    mrca(alignment_file, oldest_records_filename, copynumber_file,
         outfile, oldest_timepoint.id)


def name_to_date(name):
    id_to_date = {t.id : t.date for t in globals_.timepoints}
    return id_to_date[name.split("_")[0]]


@must_work()
def make_sequences_json(infiles, outfile):
    alignment_file, mrca_file, coords_file = infiles
    result = {}

    # add observed sequences
    for r in SeqIO.parse(alignment_file, "fasta"):
        date = name_to_date(r.name)
        if date not in result:
            result[date] = {"Observed": {}}
        result[date]["Observed"][r.id] = str(r.seq)

    # add MRCA
    mrca = read_single_record(mrca_file, 'fasta', True)
    result['MRCA'] = str(mrca.seq)

    # add reference
    reference = read_single_record(globals_.config.get('Parameters', 'reference_sequence'), 'fasta', True)
    rstr = str(reference.seq)
    with open(coords_file) as handle:
        alignment_coords = json.load(handle)['coordinates']
    ref_coords = open(globals_.config.get('Parameters', 'reference_coordinates')).read().strip().split()
    cmap = dict((c, i) for i, c in enumerate(ref_coords))
    ref_result = ''.join(list(rstr[cmap[c]] for c in alignment_coords))
    result['Reference'] = ref_result

    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


@must_work()
def make_trees_json(infile, outfile):
    with open(infile) as handle:
        newick_string = handle.read()
    result = {'Combined':
                  {'all':
                       {'Maximum Likelihood': newick_string}}}
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


def extend_coordinates(coordinates, seq, gap=None):
    """Extend coordinates to a gappy sequence.

    >>> transfer_coordinates([1, 2, 3, 4], "a-b-cd")
    [1, 1, 2, 2, 3, 4]

    """
    if gap is None:
        gap = "-"
    if sum(1 for c in seq if c != gap) != len(coordinates):
        raise Exception('coordinates do not match source')
    coords_gen = iter(coordinates)
    cur = -1
    result = []
    for char in seq:
        if char != gap:
            cur = next(coords_gen)
        result.append(cur)
    return result


@must_work()
def make_coordinates_json(infile, outfile):
    # FIXME: degap MRCA before running?
    # FIXME: split up so mafft can run remotely
    combined = os.path.join(pipeline_dir, 'combined.fasta')
    aligned = os.path.join(pipeline_dir, 'combined.aligned.fasta')
    cat_files([infile, globals_.config.get('Parameters', 'reference_sequence')], combined)
    mafft(combined, aligned)
    pairwise_mrca, pairwise_ref = list(SeqIO.parse(aligned, 'fasta'))
    ref_coords = open(globals_.config.get('Parameters', 'reference_coordinates')).read().strip().split()

    # transfer coordinates to MRCA
    pairwise_coords = extend_coordinates(ref_coords, str(pairwise_ref.seq))
    mrca_coords = list(coord for char, coord in zip(str(pairwise_mrca.seq),
                                                    pairwise_coords)
                       if char != "-")

    # extend MRCA coords to full MSA
    mrca = read_single_record(infile, 'fasta', True)
    result = extend_coordinates(mrca_coords, str(mrca.seq))

    rdict = {'coordinates': result}
    with open(outfile, 'w') as handle:
        json.dump(rdict, handle, separators=(",\n", ":"))


@must_work()
# TODO: modify turnover script to not need json input
def make_frequencies_json(infile, outfile):
    id_to_date = {t.id : t.date for t in globals_.timepoints}
    def name_to_date(name):
        return id_to_date[name.split("_")[0]]

    records = SeqIO.parse(infile, "fasta")
    result = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for r in records:
        for i, residue in enumerate(str(r.seq)):
            result[i + 1][name_to_date(r.id)][residue] += 1
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


@must_work()
def turnover(infile, outfile):
    script_path = os.path.join(globals_.script_dir, 'AAturn.jl')
    call("julia {script} {infile} {outfile} 0.01".format(script=script_path,
                                                         infile=infile, outfile=outfile))


@must_work()
def make_rates_json(infile, outfile):
    # FIXME: this is just a dummy file.
    first = read_single_record(infile, 'fasta')
    n_posns = len(first)
    dummy_rates = [[0, 0, 0, 0]] * n_posns
    result = {t.date : dummy_rates for t in globals_.timepoints}
    result['Combined'] = dummy_rates
    with open(outfile, 'w') as handle:
        json.dump(result, handle)


def make_fasttree_pipeline(name=None):
    """Factory for the FastTree sub-pipeline."""
    if name is None:
        name = "fasttree_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    translate_task = pipeline.merge(gapped_translate_wrapper,
                                    name='translate_gapped',
                                    input=None,
                                    output=os.path.join(pipeline_dir, 'translated.fasta'))
    translate_task.jobs_limit(n_local_jobs, local_job_limiter)
    # all head tasks need to make this directory, because their run order is unspecified
    translate_task.mkdir(pipeline_dir)

    fasttree_task = pipeline.merge(run_fasttree,
                                   input=None,
                                   output=os.path.join(pipeline_dir, 'tree.newick'))
    fasttree_task.jobs_limit(n_remote_jobs, remote_job_limiter)
    fasttree_task.mkdir(pipeline_dir)

    mrca_task = pipeline.merge(compute_mrca,
                               name='fasttree_mrca',
                               input=None,
                               output=os.path.join(pipeline_dir, 'mrca.fasta'))
    mrca_task.jobs_limit(n_local_jobs, local_job_limiter)
    mrca_task.mkdir(pipeline_dir)

    translate_mrca_task = pipeline.transform(gapped_translate_wrapper,
                                             name='translate_mrca',
                                             input=mrca_task,
                                             filter=suffix('.fasta'),
                                             output='.translated.fasta')
    translate_mrca_task.jobs_limit(n_local_jobs, local_job_limiter)

    coordinates_json_task = pipeline.transform(make_coordinates_json,
                                               input=translate_mrca_task,
                                               filter=formatter(),
                                               output=os.path.join(pipeline_dir, 'coordinates.json'))
    coordinates_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    sequences_json_task = pipeline.merge(make_sequences_json,
                                         input=[translate_task, translate_mrca_task, make_coordinates_json],
                                         output=os.path.join(pipeline_dir, 'sequences.json'))
    sequences_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    trees_json_task = pipeline.transform(make_trees_json,
                                         input=fasttree_task,
                                         filter=formatter(),
                                         output=os.path.join(pipeline_dir, 'trees.json'))
    trees_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    frequencies_task = pipeline.transform(make_frequencies_json,
                                               input=translate_task,
                                               filter=formatter(),
                                               output=os.path.join(pipeline_dir, 'frequencies.json'))
    frequencies_task.jobs_limit(n_local_jobs, local_job_limiter)

    turnover_task = pipeline.transform(turnover,
                                       input=frequencies_task,
                                       filter=formatter(),
                                       output=os.path.join(pipeline_dir, 'turnover.json'))
    turnover_task.jobs_limit(n_local_jobs, local_job_limiter)

    rates_json_task = pipeline.transform(make_rates_json,
                                         input=translate_task,
                                         filter=formatter(),
                                         output=os.path.join(pipeline_dir, 'rates.json'))
    rates_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    pipeline.set_head_tasks([translate_task, fasttree_task, mrca_task])
    return pipeline
