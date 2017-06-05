"""Prepare input as though it came from the consensus pipeline.

This may go before either the alignment or the analysis pipelines.

"""
import os
from shutil import copyfile

from ruffus import Pipeline, formatter

from Bio import SeqIO

import flea_pipeline.pipeline_globals as globals_
from flea_pipeline.util import must_work
from flea_pipeline.util import report_wrapper
from flea_pipeline.util import split_name


pipeline_dir = os.path.join(globals_.results_dir, "preanalysis")


def rename_record(record):
    key, rest = split_name(record.id)
    label = globals_.key_to_label[key]
    new_name = "{}_{}".format(label, rest.strip("_"))
    record = record[:]
    record.id = new_name
    record.name = ""
    record.description = ""
    return record


@must_work()
@report_wrapper
def rename_records(infile, outfile):
    records = SeqIO.parse(infile, 'fasta')
    result = (rename_record(r) for r in records)
    SeqIO.write(result, outfile, 'fasta')


@must_work()
@report_wrapper
def make_copynumbers(infile, outfile):
    cn_file = globals_.options.copynumbers
    if cn_file is None:
        records = SeqIO.parse(infile, 'fasta')
        with open(outfile, 'w') as handle:
            for r in records:
                handle.write('{}\t{}\n'.format(r.id, 1))
    else:
        copyfile(cn_file, outfile)


def make_preanalysis_pipeline(name=None):
    if name is None:
        name = "preanalysis_pipeline"
    pipeline = Pipeline(name)

    rename_records_task = pipeline.transform(rename_records,
                                             name="rename_records",
                                             input=None,
                                             filter=formatter(),
                                             output=os.path.join(pipeline_dir, "input_renamed.fasta"))

    make_copynumbers_task = pipeline.transform(make_copynumbers,
                                               name="make_copynumbers",
                                               input=rename_records,
                                               filter=formatter(),
                                               output=os.path.join(pipeline_dir, "copynumbers.tsv"))

    pipeline.set_head_tasks([rename_records_task])
    pipeline.set_tail_tasks([rename_records_task,
                             make_copynumbers_task])

    for task in pipeline.head_tasks:
        task.mkdir(pipeline_dir)
    return pipeline
