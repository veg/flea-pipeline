import os
from shutil import copyfile

from ruffus import Pipeline, suffix

from Bio import SeqIO

import pipeline_globals as globals_
from util import must_work, n_jobs, local_job_limiter, remote_job_limiter
from util import split_name


pipeline_dir = os.path.join(globals_.data_dir, "preanalysis")


def rename_record(record):
    key, rest = split_name(record.id)
    label = globals_.key_to_label[key]
    new_name = "{}_{}".format(label, rest.strip("_"))
    record.name = new_name
    record.id = new_name
    return record


@must_work()
def rename_records(infile, outfile):
    records = SeqIO.parse(infile, 'fasta')
    result = (rename_record(r) for r in records)
    SeqIO.write(result, outfile, 'fasta')


@must_work()
def make_copynumbers(infiles, outfile):
    alignfile, cnfile = infiles
    if cnfile is None:
        records = SeqIO.parse(alignfile, 'fasta')
        with open(outfile, 'w') as handle:
            for r in records:
                handle.write('{}\t{}\n'.format(r.id, 1))
    else:
        copyfile(cnfile, outfile)


def make_preanalysis_pipeline(copynumbers_file, name=None):
    if name is None:
        name = "preanalysis_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    rename_records_task = pipeline.merge(rename_records,
                                         name="rename_records",
                                         input=None,
                                         output=os.path.join(pipeline_dir, "msa_renamed.fasta"))

    make_copynumbers_task = pipeline.merge(make_copynumbers,
                                           name="make_copynumbers",
                                           input=[rename_records, copynumbers_file],
                                           output=os.path.join(pipeline_dir, "copynumbers.tsv"))
    make_copynumbers_task.jobs_limit(n_local_jobs, local_job_limiter)

    pipeline.set_head_tasks([rename_records_task])
    pipeline.set_tail_tasks([rename_records_task,
                             make_copynumbers_task])

    for task in pipeline.head_tasks:
        task.mkdir(pipeline_dir)
    return pipeline
