import os
from collections import defaultdict
import warnings

from ruffus import Pipeline, suffix, formatter, add_inputs
from Bio import SeqIO

from flea_pipeline.util import run_command
from flea_pipeline.util import new_record_seq_str, insert_gaps, update_record_seq
from flea_pipeline.util import must_work, report_wrapper
from flea_pipeline.util import check_suffix
from flea_pipeline.util import grouper
from flea_pipeline.util import translate_helper
from flea_pipeline.util import usearch_hqcs_ids

import flea_pipeline.pipeline_globals as globals_

from flea_pipeline.consensus_pipeline import cat_wrapper_ids


pipeline_dir = os.path.join(globals_.data_dir, "diagnosis")


@must_work()
@report_wrapper
def make_inputs(infiles, outfiles):
    if not len(infiles) == len(outfiles):
        raise Exception('infiles did not match outfiles')
    # ensure ccs files match. assumes first part of the name is the label.
    # FIXME: code duplication
    infiles[3:] = list(sorted(infiles[3:]))
    outfiles[3:] = list(sorted(outfiles[3:]))
    for i, o in zip(infiles, outfiles):
        if os.path.exists(o):
            os.unlink(o)
        os.symlink(i, o)


@must_work(illegal_chars='-')
@report_wrapper
def degap(infile, outfile):
    records = (SeqIO.parse(infile, 'fasta'))
    processed = (update_record_seq(r, r.seq.ungap('-')) for r in records)
    SeqIO.write(processed, outfile, 'fasta')


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
    check_suffix(infile, '.ccs.fastq')
    check_suffix(dbfile, '.degapped.fasta')
    identity = globals_.config.get('Parameters', 'ccs_to_hqcs_identity')
    maxqt = globals_.config.get('Parameters', 'max_ccs_length_ratio')
    usearch_hqcs_ids(infile, outfile, dbfile, identity,
                     maxqt, name='hqcs_ccs_ids')


@must_work()
@report_wrapper
def combine_pairs(infiles, outfiles, outdir):
    # FIXME: code duplication with compute_copynumbers
    for f in outfiles:
        os.unlink(f)
    ccsfile, hqcsfile, pairfile  = infiles
    check_suffix(ccsfile, '.ccs.fastq')
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
    stdout = os.path.join(globals_.job_script_dir, '{}.stdout'.format(outfile))
    stderr = os.path.join(globals_.job_script_dir, '{}.stderr'.format(outfile))
    return run_command(cmd, infile, outfiles=outfile, stdout=stdout, stderr=stderr)


@must_work()
@report_wrapper
def convert_bam_to_fasta(infile, outfile):
    cmd = "{} {} {}".format(globals_.config.get('Paths', 'bam2msa'), infile, outfile)
    stdout = os.path.join(globals_.job_script_dir, '{}.stdout'.format(outfile))
    stderr = os.path.join(globals_.job_script_dir, '{}.stderr'.format(outfile))
    return run_command(cmd, infile, outfiles=outfile, stdout=stdout, stderr=stderr)


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
    infile, backtranslated = infiles
    check_suffix(infile, '.aligned.fasta')
    check_suffix(backtranslated, 'codon_alignment.fasta')
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
    ccs, hqcs, cn = infiles
    check_suffix(ccs, '.no-partial-gaps.translated.fasta')
    check_suffix(hqcs, 'protein_alignment.fasta')
    check_suffix(cn, 'copynumbers.tsv')
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "diagnose.py"),
        'hqcs': hqcs,
        'ccs': ccs,
        'cn': cn,
        'd': os.path.join(pipeline_dir, 'results'),
        }
    cmd = "{python} {script} {hqcs} {ccs} {cn} {d}".format(**kwargs)
    stdout = os.path.join(globals_.job_script_dir, 'diagnose.stdout')
    stderr = os.path.join(globals_.job_script_dir, 'diagnose.stderr')
    return run_command(cmd, infiles, outfiles=outfiles,
                      stdout=stdout, stderr=stderr,
                      name="diagnose-alignment")


def make_diagnosis_pipeline(name=None):
    if name is None:
        name = "diagnosis_pipeline"
    pipeline = Pipeline(name)

    innames = [
            'protein_alignment.fasta',
            'codon_alignment.fasta',
            'copynumbers.tsv',
            ]
    inputs = list(os.path.join(pipeline_dir, f) for f in innames)
    indict = dict((k, v) for k, v in zip(innames, inputs))

    indict['ccs_sequences'] = list(os.path.join(pipeline_dir, '{}.ccs.fastq'.format(t.label))
                                   for t in globals_.timepoints)
    inputs.extend(indict['ccs_sequences'])

    make_inputs_task = pipeline.merge(make_inputs,
                                      input=None,
                                      output=inputs)
    make_inputs_task.mkdir(pipeline_dir)

    degap_backtranslated_alignment_task = pipeline.transform(degap,
                                                             name='degap_backtranslated',
                                                             input=indict['codon_alignment.fasta'],
                                                             filter=formatter('.*/(?P<LABEL>.+).fasta'),
                                                             output=os.path.join(pipeline_dir, '{LABEL[0]}.degapped.fasta'))
    degap_backtranslated_alignment_task.follows(make_inputs_task)

    hqcs_ccs_pairs_task = pipeline.transform(hqcs_ccs_pairs,
                                             input=indict['ccs_sequences'],
                                             add_inputs=add_inputs(degap_backtranslated_alignment_task),
                                             filter=formatter('.*/(?P<LABEL>.+).ccs'),
                                             output=os.path.join(pipeline_dir, '{LABEL[0]}.hqcs-ccs-pairs.txt'))

    combine_pairs_task = pipeline.subdivide(combine_pairs,
                                            input=indict['ccs_sequences'],
                                            filter=formatter('.*/(?P<LABEL>.+).ccs'),
                                            add_inputs=add_inputs(degap_backtranslated_alignment_task,
                                                                  os.path.join(pipeline_dir, '{LABEL[0]}.hqcs-ccs-pairs.txt')),
                                            output=os.path.join(pipeline_dir,
                                                                '{LABEL[0]}.ccs-alignments/combined.*.unaligned.fasta'),
                                            extras=[os.path.join(pipeline_dir, '{LABEL[0]}.ccs-alignments')])
    combine_pairs_task.mkdir(hqcs_ccs_pairs_task,
                             formatter('.*/(?P<LABEL>.+).hqcs-ccs-pairs.txt'),
                             '{path[0]}/{LABEL[0]}.ccs-alignments')

    codon_align_task = pipeline.transform(codon_align,
                                          input=combine_pairs_task,
                                          filter=suffix('.unaligned.fasta'),
                                          output='.aligned.bam')

    convert_bam_to_fasta_task = pipeline.transform(convert_bam_to_fasta,
                                                   input=codon_align_task,
                                                   filter=suffix('.bam'),
                                                   output='.fasta')

    insert_gaps_task = pipeline.transform(insert_gaps_wrapper,
                                          name='insert_gaps',
                                          input=convert_bam_to_fasta_task,
                                          add_inputs=add_inputs(indict['codon_alignment.fasta']),
                                          filter=formatter('.*/(?P<LABEL>.+).fasta'),
                                          output='{path[0]}/{LABEL[0]}.gapped.fasta')

    merge_all_timepoints_task = pipeline.merge(cat_wrapper_ids,
                                               name='merge_all_timepoints',
                                               input=insert_gaps_task,
                                               output=os.path.join(pipeline_dir, 'ccs.aligned.fasta'))

    replace_gapped_codons_task = pipeline.transform(replace_gapped_codons_file,
                                                    name='replace_gapped_codons',
                                                    input=merge_all_timepoints_task,
                                                    filter=suffix('.fasta'),
                                                    output='.no-partial-gaps.fasta')

    translate_ccs_task = pipeline.transform(gapped_translate_wrapper,
                                            name='translate_ccs',
                                            input=replace_gapped_codons_task,
                                            filter=suffix('.fasta'),
                                            output='.translated.fasta')

    diagnosis_output = list(os.path.join(pipeline_dir,
                                         'results',
                                         "freq_agreement_no_x_{}.png".format(t.label))
                            for t in globals_.timepoints)
    diagnose_alignment_task = pipeline.transform(diagnose_alignment,
                                                 name='diagnose_alignment',
                                                 input=translate_ccs_task,
                                                 add_inputs=add_inputs(indict['protein_alignment.fasta'],
                                                                       indict['copynumbers.tsv']),
                                                 filter=formatter(),
                                                 output=diagnosis_output)
    diagnose_alignment_task.mkdir(os.path.join(pipeline_dir, 'results'))

    pipeline.set_head_tasks([make_inputs_task])

    return pipeline
