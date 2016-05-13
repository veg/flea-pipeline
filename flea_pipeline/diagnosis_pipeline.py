import os
import re
from random import sample
import tempfile
from functools import partial
from collections import defaultdict
import itertools
import warnings
import shutil

import numpy as np
from ruffus import Pipeline, suffix, formatter, add_inputs
from Bio import SeqIO

from flea_pipeline.util import maybe_qsub, cat_files, touch, strlist, traverse
from flea_pipeline.util import new_record_seq_str, insert_gaps, update_record_seq
from flea_pipeline.util import must_work, must_produce, report_wrapper
from flea_pipeline.util import local_job_limiter, remote_job_limiter
from flea_pipeline.util import check_suffix, check_basename, n_jobs
from flea_pipeline.util import read_single_record
from flea_pipeline.util import partition
from flea_pipeline.util import grouper
from flea_pipeline.util import run_regexp
from flea_pipeline.util import translate_helper
from flea_pipeline.util import remove_suffix
from flea_pipeline.util import usearch_hqcs_ids

import flea_pipeline.pipeline_globals as globals_

from flea_pipeline.consensus_pipeline import cat_wrapper_ids


pipeline_dir = os.path.join(globals_.data_dir, "diagnosis")


@must_work(illegal_chars='-')
@report_wrapper
def degap(infile, outfile):
    records = (SeqIO.parse(infile, 'fasta'))
    processed = (update_record_seq(r, r.seq.ungap('-')) for r in records)
    SeqIO.write(processed, outfile, 'fasta')


@must_work()
@report_wrapper
def make_hqcs_db(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    return maybe_qsub(cmd, infile, outfiles=outfile, name='make_hqcs_db')


# making wrappers like this is necessary because nested function
# definitions are not picklable.

@must_work()
@report_wrapper
def gapped_translate_wrapper(infile, outfile):
    return translate_helper(infile, outfile, gapped=True, name='translate-gapped')


@must_work()
@report_wrapper
def hqcs_ccs_pairs(infiles, outfile):
    # FIXME: this does basically the same thing as the copynumber task,
    # except it allows CCSs to map to HQCSs in different timepoints
    infile, dbfile = infiles
    check_suffix(infile, '.lenfilter.fastq')
    check_suffix(dbfile, '.degapped.udb')
    usearch_hqcs_ids(infile, outfile, dbfile, name='hqcs_ccs_ids')


@must_work()
@report_wrapper
def combine_pairs(infiles, outfiles, outdir):
    # FIXME: code duplication with compute_copynumbers
    for f in outfiles:
        os.unlink(f)
    ccsfile, hqcsfile, pairfile = infiles
    check_suffix(ccsfile, '.lenfilter.fastq')
    check_suffix(hqcsfile, '.degapped.fasta')
    check_suffix(pairfile, '.hqcs-ccs-pairs.txt')

    with open(pairfile) as handle:
        pairs = list(line.strip().split("\t") for line in handle.readlines())
    match_dict = defaultdict(list)
    for pair in pairs:
        if len(pair) != 2:
            warnings.warn('CCS {} did not match any HQCS'.format(pair))
            continue
        ccs, hqcs = pair
        match_dict[hqcs].append(ccs)
    hqcs_records = list(SeqIO.parse(hqcsfile, "fasta"))
    ccs_records = list(SeqIO.parse(ccsfile, "fastq"))
    hqcs_dict = {r.id : r for r in hqcs_records}
    ccs_dict = {r.id : r for r in ccs_records}

    for hqcs_id, ccs_ids in match_dict.items():
        outfile = os.path.join(outdir, 'combined.{}.unaligned.fasta'.format(hqcs_id))
        records = [hqcs_dict[hqcs_id]]
        records.extend(list(ccs_dict[i] for i in ccs_ids))
        SeqIO.write(records, outfile, "fasta")


@must_work()
@report_wrapper
def codon_align(infile, outfile):
    cmd = "{} -R -a codon {} {}".format(
        globals_.config.get('Paths', 'bealign'), infile, outfile)
    stdout = os.path.join(globals_.qsub_dir, '{}.stdout'.format(outfile))
    stderr = os.path.join(globals_.qsub_dir, '{}.stderr'.format(outfile))
    return maybe_qsub(cmd, infile, outfiles=outfile, stdout=stdout, stderr=stderr)


@must_work()
@report_wrapper
def convert_bam_to_fasta(infile, outfile):
    cmd = "{} {} {}".format(globals_.config.get('Paths', 'bam2msa'), infile, outfile)
    stdout = os.path.join(globals_.qsub_dir, '{}.stdout'.format(outfile))
    stderr = os.path.join(globals_.qsub_dir, '{}.stderr'.format(outfile))
    return maybe_qsub(cmd, infile, outfiles=outfile, stdout=stdout, stderr=stderr)


def handle_gap_codon(codon):
    codon = ''.join(codon)
    if '-' in codon and codon != '---':
        return 'NNN'
    return codon


def replace_gapped_codons(record):
    if len(record.seq) % 3 != 0:
        raise Exception('record {} is not in frame'.format(record.id))
    new_str = "".join(handle_gap_codon(codon) for codon in grouper(record.seq, 3))
    result = new_record_seq_str(record, new_str)
    # HyPhy does not like anything but the sequence id
    result.name = ""
    result.description = ""
    return result


@must_work()
@report_wrapper
def replace_gapped_codons_file(infile, outfile):
    records = list(SeqIO.parse(infile, 'fasta'))
    result = (replace_gapped_codons(record) for record in records)
    SeqIO.write(result, outfile, 'fasta')


@must_work()
@report_wrapper
def insert_gaps_wrapper(infiles, outfile):
    backtranslated, infile = infiles
    check_suffix(infile, '.aligned.fasta')
    check_suffix(backtranslated, '.backtranslated.fasta')
    ref, *seqs = list(SeqIO.parse(infile, 'fasta'))
    ref_gapped = next(r for r in SeqIO.parse(backtranslated, 'fasta')
                      if r.id == ref.id)
    seqs_gapped = (new_record_seq_str(r, insert_gaps(str(ref_gapped.seq),
                                                     str(r.seq),
                                                     '-', '-'))
                   for r in seqs)
    SeqIO.write(seqs_gapped, outfile, "fasta")


@must_work()
@report_wrapper
def diagnose_alignment(infiles, outfiles):
    hqcs, cn, ccs = infiles
    check_suffix(hqcs, 'translated.aligned.edited.fasta')
    check_suffix(ccs, '.no-partial-gaps.translated.fasta')
    check_suffix(cn, '.tsv')
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "diagnose.py"),
        'hqcs': hqcs,
        'ccs': ccs,
        'cn': cn,
        'd': os.path.join(pipeline_dir, 'diagnosis')
        }
    cmd = "{python} {script} {hqcs} {ccs} {cn} {d}".format(**kwargs)
    stdout = os.path.join(globals_.qsub_dir, 'diagnose.stdout')
    stderr = os.path.join(globals_.qsub_dir, 'diagnose.stderr')
    return maybe_qsub(cmd, infiles, outfiles=outfiles,
                      stdout=stdout, stderr=stderr,
                      name="diagnose-alignment")


def make_diagnosis_pipeline(name=None):
    if name is None:
        name = "diagnosis_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    # expects backtranslated alignment as input
    degap_backtranslated_alignment_task = pipeline.transform(degap,
                                                             name='degap_backtranslated',
                                                             input=None,
                                                             filter=formatter('.*/(?P<LABEL>.+).fasta'),
                                                             output=os.path.join(pipeline_dir, '{LABEL[0]}.degapped.fasta'))
    degap_backtranslated_alignment_task.jobs_limit(n_local_jobs, local_job_limiter)
    degap_backtranslated_alignment_task.mkdir(pipeline_dir)

    make_hqcs_db_task = pipeline.transform(make_hqcs_db,
                                           input=degap_backtranslated_alignment_task,
                                           filter=suffix('.fasta'),
                                           output='.udb')
    make_hqcs_db_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    # expects CCS sequences as input
    hqcs_ccs_pairs_task = pipeline.transform(hqcs_ccs_pairs,
                                             input=None,
                                             add_inputs=add_inputs(make_hqcs_db),
                                             filter=formatter('.*/(?P<LABEL>.+).qfilter'),
                                             output=os.path.join(pipeline_dir, '{LABEL[0]}.hqcs-ccs-pairs.txt'))
    hqcs_ccs_pairs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    # expects CCS sequences as input
    combine_pairs_task = pipeline.subdivide(combine_pairs,
                                            input=None,
                                            filter=formatter('.*/(?P<LABEL>.+).qfilter'),
                                            add_inputs=add_inputs(degap_backtranslated_alignment_task,
                                                                  os.path.join(pipeline_dir, '{LABEL[0]}.hqcs-ccs-pairs.txt')),
                                            output=os.path.join(pipeline_dir,
                                                                '{LABEL[0]}.ccs-alignments/combined.*.unaligned.fasta'),
                                            extras=[os.path.join(pipeline_dir, '{LABEL[0]}.ccs-alignments')])
    combine_pairs_task.jobs_limit(n_local_jobs, local_job_limiter)
    combine_pairs_task.mkdir(hqcs_ccs_pairs_task,
                             formatter('.*/(?P<LABEL>.+).hqcs-ccs-pairs.txt'),
                             '{path[0]}/{LABEL[0]}.ccs-alignments')

    codon_align_task = pipeline.transform(codon_align,
                                          input=combine_pairs_task,
                                          filter=suffix('.unaligned.fasta'),
                                          output='.aligned.bam')
    codon_align_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    convert_bam_to_fasta_task = pipeline.transform(convert_bam_to_fasta,
                                                   input=codon_align_task,
                                                   filter=suffix('.bam'),
                                                   output='.fasta')
    convert_bam_to_fasta_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    # expects backtranslated alignment in add_inputs()
    insert_gaps_task = pipeline.transform(insert_gaps_wrapper,
                                          name='insert_gaps',
                                          input=None,
                                          add_inputs=add_inputs(convert_bam_to_fasta_task),
                                          filter=formatter('.*/(?P<LABEL>.+).fasta'),
                                          output=os.path.join(pipeline_dir, '{LABEL[0]}.gapped.fasta'))

    insert_gaps_task.jobs_limit(n_local_jobs, local_job_limiter)

    merge_all_timepoints_task = pipeline.merge(cat_wrapper_ids,
                                               name='merge_all_timepoints',
                                               input=insert_gaps_task,
                                               output=os.path.join(pipeline_dir, 'ccs.aligned.fasta'))
    merge_all_timepoints_task.jobs_limit(n_local_jobs, local_job_limiter)

    replace_gapped_codons_task = pipeline.transform(replace_gapped_codons_file,
                                                    name='replace_gapped_codons',
                                                    input=merge_all_timepoints_task,
                                                    filter=suffix('.fasta'),
                                                    output='.no-partial-gaps.fasta')
    replace_gapped_codons_task.jobs_limit(n_local_jobs, local_job_limiter)

    translate_ccs_task = pipeline.transform(gapped_translate_wrapper,
                                            name='translate_ccs',
                                            input=replace_gapped_codons_task,
                                            filter=suffix('.fasta'),
                                            output='.translated.fasta')
    translate_ccs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    # expects protein alignment and copynumbers as added input
    diagnosis_output = list(os.path.join(pipeline_dir,
                                         'diagnosis',
                                         "freq_agreement_no_x_{}.png".format(t.label))
                            for t in globals_.timepoints)
    diagnose_alignment_task = pipeline.transform(diagnose_alignment,
                                                 name='diagnose_alignment',
                                                 input=None,
                                                 add_inputs=add_inputs(translate_ccs_task),
                                                 filter=formatter(),
                                                 output=diagnosis_output)
    diagnose_alignment_task.jobs_limit(n_remote_jobs, remote_job_limiter)
    diagnose_alignment_task.mkdir(os.path.join(pipeline_dir, 'diagnosis'))

    return pipeline
