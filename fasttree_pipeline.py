"""HyPhy-free analysis. Currently only generates trees and sequences."""

import os
import json

from ruffus import Pipeline, formatter, suffix

from translate import translate
import pipeline_globals as globals_
from util import must_work, maybe_qsub, call, n_jobs, local_job_limiter, remote_job_limiter


pipeline_dir = os.path.join(globals_.data_dir, "fasttree")


@must_work()
def gapped_translate_wrapper(infile, outfile):
    translate(infile, outfile, gapped=True)


@must_work()
def run_fasttree(infile, outfile):
    stderr = '{}.stderr'.format(outfile)
    cmd = 'FastTree {} > {} 2>{}'.format(infile, outfile, stderr)
    maybe_qsub(cmd, outfiles=outfile, stdout='/dev/null', stderr='/dev/null')


def make_fasttree_pipeline(name=None):
    """Factory for the FastTree sub-pipeline."""
    if name is None:
        name = "fasttree_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    translate_task = pipeline.transform(gapped_translate_wrapper,
                                        name='translate_for_fasttree',
                                        input=None,
                                        filter=formatter(),
                                        output=os.path.join(pipeline_dir, 'translated.fasta'))
    translate_task.jobs_limit(n_local_jobs, local_job_limiter)
    translate_task.mkdir(pipeline_dir)

    fasttree_task = pipeline.transform(run_fasttree,
                                       input=translate_task,
                                       filter=suffix('.fasta'),
                                       output='.tree')
    fasttree_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    pipeline.set_head_tasks([translate_task])
    return pipeline
