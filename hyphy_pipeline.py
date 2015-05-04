import os
import json
from datetime import datetime
import shutil

from ruffus import Pipeline, formatter

from Bio import SeqIO

import pipeline_globals as globals_
from DNAcons import consfile
from util import must_work, maybe_qsub, call, n_jobs, local_job_limiter, remote_job_limiter


hyphy_script_dir = os.path.join(globals_.script_dir, 'hyphy_scripts')
hyphy_data_dir = os.path.join(globals_.data_dir, "hyphy_data")
hyphy_input_dir = os.path.join(hyphy_data_dir, "input")
hyphy_results_dir = os.path.join(hyphy_data_dir, "results")


def hyphy_script(name):
    return os.path.join(hyphy_script_dir, name)


def hyphy_input(s):
    return os.path.join(hyphy_input_dir, s)


def hyphy_results(s):
    return os.path.join(hyphy_results_dir, s)


@must_work()
def copy_rename_alignment(infile, outfile):
    # rename sequence ids to match those expected by HYPHY:
    # VXX_number_copynumber
    records = list(SeqIO.parse(infile, 'fasta'))
    for r in records:
        found = False
        # check every possible key, in case they are different lengths.
        for k, v in globals_.timepoint_ids.items():
            if r.id.startswith(k):
                rest = r.id[len(k):]
                rest = rest.strip('_')
                if not rest:
                    raise Exception('sequence id has no unique part')
                r.id = '{}_{}_1'.format(v, rest)
                r.name = ''
                r.description = ''
                found = True
                break
        if not found:
            raise Exception('record id not found in input file')
    SeqIO.write(records, outfile, 'fasta')


@must_work()
def copy_wrapper(infile, outfile):
    shutil.copyfile(infile, outfile)


@must_work()
def write_dates(infile, outfile):
    # NOTE: we assume the first part of the record id is the timestamp
    # id, followed by an underscore.
    records = list(SeqIO.parse(infile, "fasta"))
    id_to_date = {t.id : t.date for t in globals_.timepoints}
    with open(outfile, "w") as handle:
        outdict = {r.id: id_to_date[r.id.split("_")[0]] for r in records}
        json.dump(outdict, handle, separators=(",\n", ":"))


def mrca(infile, recordfile, outfile, oldest_id):
    """Writes records from `infile` to `recordfile` that have an ID
    corresponding to `oldest_id`. Then runs DNAcons, writing result to
    `outfile`.

    """
    len_id = len(oldest_id)
    records = SeqIO.parse(infile, "fasta")
    oldest_records = (r for r in records if r.id.startswith(oldest_id))
    SeqIO.write(oldest_records, recordfile, "fasta")
    consfile(recordfile, outfile, id_str="mrca", ungap=False)


@must_work(illegal_chars='')
def compute_mrca(infile, outfile):
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(globals_.timepoints, key=strptime)
    oldest_records_filename = '.'.join([infile, "oldest_{}".format(oldest_timepoint.id)])
    mrca(infile, oldest_records_filename, outfile, oldest_timepoint.id)


def hyphy_call(script_file, name, args):
    if args:
        in_str = "".join("{}\n".format(i) for i in args)
    else:
        in_str = ''
    infile = '{}.stdin'.format(name)
    with open(infile, 'w') as handle:
        handle.write(in_str)
    cmd = '{} {} < {}'.format(globals_.config.get('Paths', 'hyphy'), script_file, infile)
    maybe_qsub(cmd, name=name)


@must_work()
def compute_hxb2_coords(infile, outfile):
    hxb2_script = hyphy_script('HXB2partsSplitter.bf')
    hyphy_call(hxb2_script, 'hxb2_coords', [infile, outfile])


@must_work()
def evo_history(infiles, outfile):
    hyphy_call(hyphy_script("obtainEvolutionaryHistory.bf"),
               'evo_history', [hyphy_data_dir])


@must_work()
def aa_freqs(infile, outfile):
    hyphy_call(hyphy_script("aminoAcidFrequencies.bf"),
               'aa_freqs',
               [hyphy_data_dir])


@must_work()
def turnover(infile, outfile):
    script_path = os.path.join(globals_.script_dir, 'AAturn.jl')
    call("julia {script} {infile} {outfile} 0.01".format(script=script_path,
                                                         infile=infile, outfile=outfile))


@must_work()
def run_fubar(infile, outfile):
    hyphy_call(hyphy_script('runFUBAR.bf'), 'fubar', [hyphy_data_dir])


def make_hyphy_pipeline(standalone, name=None):
    """Factory for the hyphy sub-pipeline."""
    if name is None:
        name = "hyphy_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    if standalone:
        copy_function = copy_rename_alignment
    else:
        copy_function = copy_wrapper

    copy_task = pipeline.transform(copy_function,
                                   name="copy_alignment",
                                   input=None,
                                   filter=formatter(),
                                   output=hyphy_input('merged.fas'))
    copy_task.jobs_limit(n_local_jobs, local_job_limiter)
    copy_task.mkdir(hyphy_input_dir)
    copy_task.mkdir(hyphy_results_dir)
    pipeline.set_head_tasks([copy_task])

    dates_task = pipeline.transform(write_dates,
                                    input=copy_task,
                                    filter=formatter(),
                                    output=hyphy_input('merged.dates'))
    dates_task.jobs_limit(n_local_jobs, local_job_limiter)

    mrca_task = pipeline.transform(compute_mrca,
                                   input=copy_task,
                                   filter=formatter(),
                                   output=hyphy_input('earlyCons.seq'))
    mrca_task.jobs_limit(n_local_jobs, local_job_limiter)

    coords_task = pipeline.transform(compute_hxb2_coords,
                                     input=mrca_task,
                                     filter=formatter(),
                                     output=hyphy_input('merged.prot.parts'))
    coords_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    evo_history_output = list(hyphy_results(f) for f in ('rates_pheno.tsv',
                                                         'trees.json',
                                                         'sequences.json'))
    evo_history_output.append(hyphy_input('mrca.seq'))
    evo_history_task = pipeline.merge(evo_history,
                                      input=[dates_task, coords_task, mrca_task, copy_task],
                                      output=evo_history_output)
    evo_history_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    aa_freqs_task = pipeline.merge(aa_freqs,
                                   input=[dates_task, mrca_task, evo_history_task, copy_task],
                                   output=hyphy_results('frequencies.json'))
    aa_freqs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    turnover_task = pipeline.transform(turnover,
                                       input=aa_freqs_task,
                                       filter=formatter(),
                                       output=hyphy_results('turnover.json'))
    turnover_task.jobs_limit(n_local_jobs, local_job_limiter)

    fubar_task = pipeline.merge(run_fubar,
                                input=[dates_task, evo_history_task, copy_task],
                                output=hyphy_results('rates.json'))
    fubar_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    return pipeline
