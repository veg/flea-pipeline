import os
import shutil

from ruffus import Pipeline, suffix, formatter, add_inputs
from Bio import SeqIO

from flea_pipeline.util import run_command, cat_files
from flea_pipeline.util import must_work, report_wrapper
from flea_pipeline.util import check_suffix, check_basename
from flea_pipeline.util import translate_helper

import flea_pipeline.pipeline_globals as globals_


pipeline_dir = os.path.join(globals_.results_dir, "alignment")


@must_work()
@report_wrapper
def make_input(infile, outfile):
    if os.path.exists(outfile):
        os.unlink(outfile)
    os.symlink(infile, outfile)


def mafft(infile, outfile):
    binary = globals_.config.get('Paths', 'mafft')
    stderr = '{}.stderr'.format(outfile)
    cmd = '{} --ep 0.5 --quiet --preservecase {} > {} 2>{}'.format(binary, infile, outfile, stderr)
    return run_command(cmd, infile, outfiles=outfile)


@must_work()
@report_wrapper
def translate_wrapper(infile, outfile):
    return translate_helper(infile, outfile, gapped=False, name='translate')


@must_work(seq_ids=True)
@report_wrapper
def mafft_wrapper_seq_ids(infile, outfile):
    mafft(infile, outfile)


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


def pause_for_editing_alignment():
    if globals_.config.getboolean('Tasks', 'pause_after_codon_alignment'):
        infile = os.path.join(pipeline_dir, 'hqcs.translated.aligned.fasta')
        outfile = os.path.join(pipeline_dir, 'hqcs.translated.aligned.edited.fasta')
        input('Paused for manual editing of {}'
              '\nPress Enter to continue.'.format(outfile))
        # check that all sequences in `outfile` differ only in gap locations
        unedited = dict((r.id, r) for r in SeqIO.parse(infile, 'fasta'))
        edited = dict((r.id, r) for r in SeqIO.parse(outfile, 'fasta'))
        for k in edited.keys():
            edited_seq = list(c for c in str(edited[k].seq) if c != '-')
            unedited_seq = list(c for c in str(unedited[k].seq) if c != '-')
            if edited_seq != unedited_seq:
                raise Exception('Incorrect manual editing of "{}".'
                                ' "{}" differs after editing.'.format(outfile, k))


@must_work()
@report_wrapper
def backtranslate_alignment(infiles, outfile):
    hqcs, aligned_protein = infiles
    check_basename(hqcs, 'hqcs.fasta')
    check_suffix(aligned_protein, '.translated.aligned.edited.fasta')
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "backtranslate.py"),
        'protein': aligned_protein,
        'dna': hqcs,
        'outfile': outfile,
        }
    cmd = "{python} {script} {protein} {dna} {outfile}".format(**kwargs)
    stdout = os.path.join(globals_.job_script_dir, 'backtranslate.stdout')
    stderr = os.path.join(globals_.job_script_dir, 'backtranslate.stderr')
    return run_command(cmd, infiles, outfiles=outfile,
                      stdout=stdout, stderr=stderr,
                      name="backtranslate")


def make_alignment_pipeline(name=None):
    if name is None:
        name = "alignment_pipeline"
    pipeline = Pipeline(name)

    inputfile = os.path.join(pipeline_dir, 'hqcs.fasta')
    make_input_task = pipeline.transform(make_input,
                                         input=None,
                                         filter=formatter(),
                                         output=inputfile)
    make_input_task.mkdir(pipeline_dir)

    translate_hqcs_task = pipeline.transform(translate_wrapper,
                                             name='translate_hqcs',
                                             input=make_input_task,
                                             filter=formatter('.*/(?P<LABEL>.+).fasta'),
                                             output=os.path.join(pipeline_dir, '{LABEL[0]}.translated.fasta'))

    align_hqcs_protein_task = pipeline.transform(mafft_wrapper_seq_ids,
                                                 name='align_hqcs_protein',
                                                 input=translate_hqcs_task,
                                                 filter=suffix('.fasta'),
                                                 output='.aligned.fasta')

    copy_protein_alignment_task = pipeline.transform(copy_file,
                                                     name='copy_protein_alignment',
                                                     input=align_hqcs_protein_task,
                                                     filter=suffix('.fasta'),
                                                     output='.edited.fasta')
    copy_protein_alignment_task.posttask(pause_for_editing_alignment)

    # expects all HQCSs as input
    backtranslate_alignment_task = pipeline.transform(backtranslate_alignment,
                                                      input=make_input_task,
                                                      add_inputs=add_inputs(copy_protein_alignment_task),
                                                      filter=formatter(),
                                                      output=os.path.join(pipeline_dir,
                                                                          'hqcs.translated.aligned.edited.backtranslated.fasta'))

    pipeline.set_head_tasks([make_input_task])
    pipeline.set_tail_tasks([backtranslate_alignment_task,
                             copy_protein_alignment_task])

    return pipeline