"""HyPhy-free analysis. Currently only generates trees and sequences."""

import os
import json
from itertools import islice

from ruffus import Pipeline, formatter, suffix

from Bio import SeqIO

from translate import translate
import pipeline_globals as globals_
from util import must_work, maybe_qsub, call, n_jobs, local_job_limiter, remote_job_limiter, read_single_record, cat_files

from alignment_pipeline import mafft
from hyphy_pipeline import compute_mrca


pipeline_dir = os.path.join(globals_.data_dir, "fasttree")


@must_work()
def gapped_translate_wrapper(infile, outfile):
    translate(infile, outfile, gapped=True)


@must_work()
def run_fasttree(infile, outfile):
    binary = globals_.config.get('Paths', 'FastTree')
    stderr = '{}.stderr'.format(outfile)
    cmd = '{} -gtr -nt {} > {} 2>{}'.format(binary, infile, outfile, stderr)
    maybe_qsub(cmd, outfiles=outfile, stdout='/dev/null', stderr='/dev/null')


def name_to_date(name):
    id_to_date = {t.id : t.date for t in globals_.timepoints}
    return id_to_date[name.split("_")[0]]


@must_work()
def make_sequences_json(infiles, outfile):
    alignment_file, mrca_file = infiles
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
# FIXME: change frontend to not require this
def make_frequencies_json(infile, outfile):
    with open(infile) as handle:
        coordinates = json.load(handle)['coordinates']
    result = {str(i + 1): {"HXB2": val} for i, val in enumerate(coordinates)}
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


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

    translate_task = pipeline.transform(gapped_translate_wrapper,
                                        name='translate_gapped',
                                        input=None,
                                        filter=formatter(),
                                        output=os.path.join(pipeline_dir, 'translated.fasta'))
    translate_task.jobs_limit(n_local_jobs, local_job_limiter)
    translate_task.mkdir(pipeline_dir)

    fasttree_task = pipeline.transform(run_fasttree,
                                       input=None,
                                       filter=formatter(),
                                       output=os.path.join(pipeline_dir, 'tree.newick'))
    fasttree_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    mrca_task = pipeline.transform(compute_mrca,
                                   name='fasttree_mrca',
                                   input=None,
                                   filter=formatter(),
                                   output=os.path.join(pipeline_dir, 'mrca.fasta'))
    mrca_task.jobs_limit(n_local_jobs, local_job_limiter)

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
                                         input=[translate_task, translate_mrca_task],
                                         output=os.path.join(pipeline_dir, 'sequences.json'))
    sequences_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    trees_json_task = pipeline.transform(make_trees_json,
                                         input=fasttree_task,
                                         filter=formatter(),
                                         output=os.path.join(pipeline_dir, 'trees.json'))
    trees_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    frequencies_json_task = pipeline.transform(make_frequencies_json,
                                               input=coordinates_json_task,
                                               filter=formatter(),
                                               output=os.path.join(pipeline_dir, 'frequencies.json'))
    frequencies_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    rates_json_task = pipeline.transform(make_rates_json,
                                         input=translate_task,
                                         filter=formatter(),
                                         output=os.path.join(pipeline_dir, 'rates.json'))
    rates_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    pipeline.set_head_tasks([translate_task, fasttree_task, mrca_task])
    return pipeline
