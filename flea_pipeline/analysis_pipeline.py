import os
import json
from itertools import islice
from datetime import datetime
from collections import defaultdict
import shutil
import csv

from ruffus import Pipeline, formatter, suffix

import numpy as np

from Bio.Seq import translate
from Bio import SeqIO
from Bio import Phylo

import flea_pipeline.pipeline_globals as globals_
from flea_pipeline.util import must_work, report_wrapper, maybe_qsub, call, n_jobs
from flea_pipeline.util import local_job_limiter, remote_job_limiter
from flea_pipeline.util import read_single_record, cat_files
from flea_pipeline.util import name_to_date
from flea_pipeline.util import extend_coordinates
from flea_pipeline.util import new_record_seq_str
from flea_pipeline.util import grouper
from flea_pipeline.util import column_count, prob, js_divergence
from flea_pipeline.util import translate_helper

from flea_pipeline.alignment_pipeline import mafft, cat_wrapper_ids


pipeline_dir = os.path.join(globals_.data_dir, "analysis")
hyphy_script_dir = globals_.script_dir


def hyphy_script(name):
    return os.path.join(hyphy_script_dir, name)


def hyphy_call(script_file, infiles, outfiles,  name, args):
    if args:
        in_str = "".join("{}\n".format(i) for i in args)
    else:
        in_str = ''
    infile = os.path.join(globals_.qsub_dir, '{}.stdin'.format(name))
    with open(infile, 'w') as handle:
        handle.write(in_str)
    cmd = '{} {} < {}'.format(globals_.config.get('Paths', 'hyphy'), script_file, infile)
    return maybe_qsub(cmd, infiles, outfiles, name=name,
                      walltime=globals_.config.getint('Jobs', 'hyphy_walltime'))


def replace_id(record, id_):
    result = record[:]
    result.id = id_
    # HyPhy fails if a name or description are present
    result.name = ""
    result.description = ""
    return result


@must_work()
@report_wrapper
def dates_json(infiles, outfile):
    result = dict((v, k) for k, v in globals_.label_to_date.items())
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


@must_work()
@report_wrapper
def copy_copynumber_file(infiles, outfile):
    _, infile = infiles
    shutil.copyfile(infile, outfile)


def parse_copynumbers(infile):
    with open(infile) as handle:
        parsed = csv.reader(handle, delimiter='\t')
        return dict((i, int(n)) for i, n in parsed)


def id_with_cn(id_, cn):
    return "{}_cn_{}".format(id_, cn)


@must_work()
@report_wrapper
def add_copynumbers(infiles, outfile):
    msafile, copyfile = infiles
    id_to_cn = parse_copynumbers(copyfile)
    records = SeqIO.parse(msafile, 'fasta')
    result = (replace_id(r, id_with_cn(r.id, id_to_cn[r.id]))
              for r in records)
    SeqIO.write(result, outfile, "fasta")


@must_work()
@report_wrapper
def copynumber_json(infile, outfile):
    d = parse_copynumbers(infile)
    # add copynumber to name, to match rest of this pipeline
    result = dict((id_with_cn(key, value), value)
                  for key, value in d.items())
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


@must_work()
@report_wrapper
def gapped_translate_wrapper(infile, outfile):
    return translate_helper(infile, outfile, gapped=True, name='translate-gapped')


@must_work()
@report_wrapper
def run_fasttree(infile, outfile):
    binary = globals_.config.get('Paths', 'FastTree')
    stderr = '{}.stderr'.format(outfile)
    cmd = '{} -gtr -nt {} > {} 2>{}'.format(binary, infile, outfile, stderr)
    return maybe_qsub(cmd, infile, outfiles=outfile, stdout='/dev/null', stderr='/dev/null')


@must_work()
@report_wrapper
def reroot_at_mrca(infile, outfile):
    tree = next(Phylo.parse(infile, 'newick'))
    clade = next(tree.find_clades('mrca'))
    tree.root_with_outgroup(clade)
    Phylo.write([tree], outfile, 'newick')


def mrca(infile, copynumber_file, outfile):
    """Runs DNAcons, writing result to `outfile`."""
    # TODO: replace with weighted alignment-free consensus
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "DNAcons.py"),
        'infile': infile,
        'outfile': outfile,
        'copynumber_file': copynumber_file,
        'ambifile': "{}.ambi".format(outfile),
        'id_str': 'mrca',
        }
    cmd = ("{python} {script} --keep-gaps --codon -o {outfile}"
           " --copynumbers {copynumber_file} --ambifile {ambifile}"
           " --id {id_str} {infile}".format(**kwargs))
    return maybe_qsub(cmd, infile, outfile, name="mrca")


@must_work(in_frame=True)
@report_wrapper
def compute_mrca(infiles, outfile):
    """Writes oldest sequences to a seperate file, then generates
    their consensus.

    """
    alignment_file, copynumber_file = infiles
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(globals_.timepoints, key=strptime)
    oldest_records_filename = os.path.join(pipeline_dir, 'oldest_sequences.fasta')

    records = SeqIO.parse(alignment_file, "fasta")
    oldest_records = (r for r in records if r.id.startswith(oldest_timepoint.label))
    SeqIO.write(oldest_records, oldest_records_filename, "fasta")

    return mrca(oldest_records_filename, copynumber_file, outfile)


@must_work()
@report_wrapper
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
@report_wrapper
def make_trees_json(infile, outfile):
    with open(infile) as handle:
        newick_string = handle.read()
    result = {'Combined':
                  {'all':
                       {'Maximum Likelihood': newick_string}}}
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


@must_work()
@report_wrapper
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
@report_wrapper
def js_divergence_json(infile, outfile):
    records = list(SeqIO.parse(infile, 'fasta'))
    seq_array = np.array(list(list(str(r.seq)) for r in records))
    # assume id has "_cn_<copynumber>" prefix
    copynumber_array = np.array(list(int(r.id.split("_cn_")[1])
                                     for r in records))
    date_array = np.array(list(name_to_date(r.id) for r in records))
    dates = list(sorted(set(date_array)))

    aas = list(sorted(set(seq_array.ravel())))
    keys = np.array(aas)

    ps = {}
    for date in dates:
        bools = (date_array == date)
        counts = column_count(seq_array[bools], keys,
                              weights=copynumber_array[bools])
        ps[date] = prob(counts, axis=0).T

    result = {}
    for i in range(1, len(dates)):
        prev_date = dates[i - 1]
        cur_date = dates[i]
        prev_ps = ps[prev_date]
        cur_ps = ps[cur_date]
        result[dates[i]] = list(js_divergence(a, b)
                                for a, b in zip(prev_ps, cur_ps))
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",", ":"))


@must_work()
@report_wrapper
def write_dates(infile, outfile):
    records = list(SeqIO.parse(infile, "fasta"))
    with open(outfile, "w") as handle:
        outdict = {r.id: name_to_date(r.id) for r in records}
        json.dump(outdict, handle, separators=(",\n", ":"))


def handle_codon(codon):
    codon = ''.join(codon)
    if codon == "---":
        return codon
    if translate(codon) == "*":
        return "NNN"
    return codon


def replace_stop_codons(record):
    if len(record.seq) % 3 != 0:
        raise Exception('record {} is not in frame'.format(record.id))
    new_str = "".join(handle_codon(codon) for codon in grouper(record.seq, 3))
    result = new_record_seq_str(record, new_str)
    # HyPhy does not like anything but the sequence id
    result.name = ""
    result.description = ""
    return result


@must_work()
@report_wrapper
def replace_stop_codons_file(infile, outfile):
    records = list(SeqIO.parse(infile, 'fasta'))
    result = (replace_stop_codons(record) for record in records)
    SeqIO.write(result, outfile, 'fasta')


@must_work()
@report_wrapper
def compute_hxb2_regions(infile, outfile):
    hxb2_script = hyphy_script('HXB2partsSplitter.bf')
    return hyphy_call(hxb2_script, infile, outfile, 'hxb2_regions', [infile, outfile])


@must_work()
@report_wrapper
def evo_history(infiles, outfile):
    outfiles = [
        os.path.join(pipeline_dir, 'unused_evo_history_mrca'),
        outfile,
        os.path.join(pipeline_dir, 'unused_evo_history_trees'),
        os.path.join(pipeline_dir, 'unused_evo_history_ancestral'),
        os.path.join(pipeline_dir, 'unused_evo_history_sequences'),
        ]
    params = infiles + outfiles
    return hyphy_call(hyphy_script("obtainEvolutionaryHistory.bf"),
                      infiles, outfile, 'evo_history', params)


@must_work()
@report_wrapper
def rates_pheno_json(infile, outfile):
    with open(infile) as h:
        lines = list(csv.reader(h, delimiter='\t'))
    result = {}
    keys, rest = lines[0], lines[1:]
    result = list(dict((key, value) for key, value in zip(keys, line))
                  for line in rest)
    with open(outfile, 'w') as h:
        json.dump(result, h, separators=(",\n", ":"))


@must_work()
@report_wrapper
def run_fubar(infiles, outfile):
    params = infiles + ["{}/".format(pipeline_dir), outfile]
    return hyphy_call(hyphy_script('runFUBAR.bf'), infiles, outfile, 'fubar', params)


def make_analysis_pipeline(do_hyphy, name=None):
    """Factory for the analysis sub-pipeline."""
    if name is None:
        name = "analysis_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    dates_json_task = pipeline.transform(dates_json,
                                         input=None,
                                         filter=formatter(),
                                         output=os.path.join(pipeline_dir, 'dates.json'))
    dates_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    # evo_history needs copynumbers in sequence names
    add_copynumbers_task = pipeline.merge(add_copynumbers,
                                          name="add_copynumbers",
                                          input=None,
                                          output=os.path.join(pipeline_dir, "msa.copynumber-ids.fasta"))
    add_copynumbers_task.jobs_limit(n_local_jobs, local_job_limiter)

    copy_copynumber_task = pipeline.merge(copy_copynumber_file,
                                          name="copy_copynumber_task",
                                          input=None,
                                          output=os.path.join(pipeline_dir, 'copynumbers.tsv'))
    copy_copynumber_task.jobs_limit(n_local_jobs, local_job_limiter)

    copynumber_json_task = pipeline.merge(copynumber_json,
                                          input=copy_copynumber_task,
                                          output=os.path.join(pipeline_dir, 'copynumbers.json'))
    copynumber_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    mrca_task = pipeline.merge(compute_mrca,
                               name='compute_mrca',
                               input=None,
                               output=os.path.join(pipeline_dir, 'mrca.fasta'))
    mrca_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    add_mrca_task = pipeline.merge(cat_wrapper_ids,
                                   name='add_mrca',
                                   input=[mrca_task, add_copynumbers_task],
                                   output=os.path.join(pipeline_dir, 'msa.with-mrca.fasta'))
    add_mrca_task.jobs_limit(n_local_jobs, local_job_limiter)

    fasttree_task = pipeline.transform(run_fasttree,
                                       input=add_mrca_task,
                                       filter=formatter(),
                                       output=os.path.join(pipeline_dir, 'tree.newick'))
    fasttree_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    reroot_task = pipeline.transform(reroot_at_mrca,
                                     input=fasttree_task,
                                     filter=suffix('.newick'),
                                     output='.rerooted.newick')

    trees_json_task = pipeline.transform(make_trees_json,
                                         input=reroot_task,
                                         filter=formatter(),
                                         output=os.path.join(pipeline_dir, 'trees.json'))
    trees_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    translate_task = pipeline.transform(gapped_translate_wrapper,
                                        name='translate_gapped',
                                        input=add_copynumbers_task,
                                        filter=suffix('.fasta'),
                                        output='.translated.fasta')
    translate_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    translate_mrca_task = pipeline.transform(gapped_translate_wrapper,
                                             name='translate_mrca',
                                             input=mrca_task,
                                             filter=suffix('.fasta'),
                                             output='.translated.fasta')
    translate_mrca_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    coordinates_json_task = pipeline.transform(make_coordinates_json,
                                               input=translate_mrca_task,
                                               filter=formatter(),
                                               output=os.path.join(pipeline_dir, 'coordinates.json'))
    coordinates_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    sequences_json_task = pipeline.merge(make_sequences_json,
                                         input=[translate_task, translate_mrca_task, make_coordinates_json],
                                         output=os.path.join(pipeline_dir, 'sequences.json'))
    sequences_json_task.jobs_limit(n_local_jobs, local_job_limiter)

    js_divergence_task = pipeline.transform(js_divergence_json,
                                            input=translate_task,
                                            filter=formatter(),
                                            output=os.path.join(pipeline_dir, 'js_divergence.json'))
    js_divergence_task.jobs_limit(n_local_jobs, local_job_limiter)

    dates_task = pipeline.transform(write_dates,
                                    input=add_copynumbers_task,
                                    filter=formatter(),
                                    output=os.path.join(pipeline_dir, 'merged.dates'))
    dates_task.jobs_limit(n_local_jobs, local_job_limiter)

    if do_hyphy:
        replace_stop_codons_task = pipeline.transform(replace_stop_codons_file,
                                                      input=add_copynumbers_task,
                                                      filter=suffix('.fasta'),
                                                      output='.no-stops.fasta')
        replace_stop_codons_task.jobs_limit(n_local_jobs, local_job_limiter)

        region_coords_task = pipeline.transform(compute_hxb2_regions,
                                                input=mrca_task,
                                                filter=formatter(),
                                                output=os.path.join(pipeline_dir, 'region_coords.json'))
        region_coords_task.jobs_limit(n_remote_jobs, remote_job_limiter)

        evo_history_task = pipeline.merge(evo_history,
                                          input=[replace_stop_codons_task, dates_task, region_coords_task, mrca_task],
                                          output=os.path.join(pipeline_dir, 'rates_pheno.tsv'))
        evo_history_task.jobs_limit(n_remote_jobs, remote_job_limiter)

        rates_pheno_json_task = pipeline.transform(rates_pheno_json,
                                                   input=evo_history_task,
                                                   filter=suffix('.tsv'),
                                                   output='.json')
        rates_pheno_json_task.jobs_limit(n_local_jobs, local_job_limiter)

        fubar_task = pipeline.merge(run_fubar,
                                    input=[replace_stop_codons_task, dates_task, mrca_task],
                                    output=os.path.join(pipeline_dir, 'rates.json'))
        fubar_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    pipeline.set_head_tasks([add_copynumbers_task, copy_copynumber_task, mrca_task, dates_json_task])
    for task in pipeline.head_tasks:
        task.mkdir(pipeline_dir)

    return pipeline
