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

import flea_pipeline.pipeline_globals as globals_


pipeline_dir = os.path.join(globals_.data_dir, "alignment")


def mafft(infile, outfile):
    stderr = '{}.stderr'.format(outfile)
    cmd = 'mafft-fftns --ep 0.5 --quiet --preservecase {} > {} 2>{}'.format(infile, outfile, stderr)
    return maybe_qsub(cmd, infile, outfiles=outfile)


def min_len():
    min_len_fraction = globals_.config.getfloat('Parameters', 'min_sequence_length')
    ref_file = globals_.config.get('Parameters', 'reference_sequence')
    # multiple by 3 because reference sequence is translated.
    return 3 * int(min_len_fraction * len(read_single_record(ref_file, 'fasta').seq))


@must_work()
@report_wrapper
def filter_fastq(infile, outfile):
    qmax = globals_.config.get('Parameters', 'qmax')
    max_err_rate = globals_.config.get('Parameters', 'max_err_rate')
    cmd = ('{usearch} -fastq_filter {infile} -fastq_maxee_rate {max_err_rate}'
         ' -threads 1 -fastq_qmax {qmax} -fastq_minlen {min_len} -fastqout {outfile}'
         ' -relabel "{seq_id}_"'.format(
            usearch=globals_.config.get('Paths', 'usearch'),
            infile=infile, outfile=outfile, qmax=qmax, min_len=min_len(),
            max_err_rate=max_err_rate, seq_id=globals_.key_to_label[infile]))
    return maybe_qsub(cmd, infile, outfiles=outfile, name='filter-fastq')


@must_work()
@report_wrapper
def trim_heads_and_tails(infile, outfile):
    python = globals_.config.get('Paths', 'python')
    script = os.path.join(globals_.script_dir, 'trim.py')
    cmd = ("{python} {script} --fastq {infile} {outfile}")
    cmd = cmd.format(python=python, script=script, infile=infile, outfile=outfile)
    return maybe_qsub(cmd, infile, outfile, name='trim-heads-and-tails')


@must_work()
@report_wrapper
def filter_runs(infile, outfile):
    runlen = globals_.config.getint('Parameters', 'run_length')
    cregexp = run_regexp(runlen)
    records = SeqIO.parse(infile, 'fastq')
    result = (r for r in records if cregexp.search(str(r.seq)) is None)
    SeqIO.write(result, outfile, 'fastq')


@must_work()
@report_wrapper
def filter_contaminants(infile, outfile):
    uncontam = outfile
    contam = outfile.replace('good', 'bad')
    outfiles = [uncontam, contam]
    cmd = ('{usearch} -usearch_global {infile} -db {db}'
           ' -id {id} -notmatched {uncontam} -matched {contam}'
           ' -strand both'.format(usearch=globals_.config.get('Paths', 'usearch'),
                                  infile=infile, db=globals_.config.get('Parameters', 'contaminants_db'),
                                  id=globals_.config.get('Parameters', 'contaminant_identity'),
                                  uncontam=uncontam, contam=contam))
    return maybe_qsub(cmd, infile, outfiles=outfiles, name='filter-contaminants')


def usearch_reference_db(infile, outfile, dbfile=None, name=None):
    """run usearch_global against reference database and print fasta hits"""
    if dbfile is None:
        dbfile = globals_.config.get('Parameters', 'reference_db')
    identity = globals_.config.get('Parameters', 'reference_identity')
    max_accepts = globals_.config.get('Parameters', 'max_accepts')
    max_rejects = globals_.config.get('Parameters', 'max_rejects')
    cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
           ' -fastapairs {outfile} -alnout {outfile}.human -userout {outfile}.calns'
           ' -userfields caln -top_hit_only -strand both'
           ' -maxaccepts {max_accepts} -maxrejects {max_rejects}')
    if name is None:
        name = 'usearch-reference-db'
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, db=dbfile, id=identity, outfile=outfile,
                     max_accepts=max_accepts, max_rejects=max_rejects)
    return maybe_qsub(cmd, infile, outfiles=outfile, name=name)


@must_work()
@report_wrapper
def filter_uncontaminated(infile, outfile):
    return usearch_reference_db(infile, outfile, name="filter-uncontaminated")


def shift_correction_helper(infile, outfile, keep, correct, name=None):
    python = globals_.config.get('Paths', 'python')
    script = os.path.join(globals_.script_dir, 'correct_shifts.py')
    keep_option = "--keep" if keep else ""
    correct_option = "--correct-dels" if correct else ""
    calns_file = "{}.calns".format(infile)
    discard_file = "{}.discarded".format(outfile)
    summary_file = "{}.summary".format(outfile)
    cmd = ("{python} {script} {keep_option} {correct_option} --calns={calns_file}"
           " --discard={discard_file} --summary={summary_file}"
           " {infile} {outfile}")
    cmd = cmd.format(python=python, script=script, keep_option=keep_option,
                     correct_option=correct_option,
                     calns_file=calns_file, discard_file=discard_file,
                     summary_file=summary_file, infile=infile, outfile=outfile)
    return maybe_qsub(cmd, infile, outfile, name=name)


@must_work(seq_ratio=(2, 1))
@report_wrapper
def ccs_shift_correction(infile, outfile):
    return shift_correction_helper(infile, outfile, keep=True, correct=False,
                                   name="ccs-shift-correction")


@must_work(seq_ratio=(2, 1))
@report_wrapper
def every_other(infile, outfile):
    records = SeqIO.parse(infile, 'fasta')
    result = itertools.islice(records, 0, None, 2)
    SeqIO.write(result, outfile, 'fasta')


@must_work(illegal_chars='-')
@report_wrapper
def degap(infile, outfile):
    records = (SeqIO.parse(infile, 'fasta'))
    processed = (update_record_seq(r, r.seq.ungap('-')) for r in records)
    SeqIO.write(processed, outfile, 'fasta')


@must_work()
@report_wrapper
def filter_length(infile, outfile):
    cutoff = min_len()
    records = SeqIO.parse(infile, 'fasta')
    result = (r for r in records if len(r.seq) >= cutoff)
    SeqIO.write(result, outfile, 'fasta')


@must_produce(n=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def cluster(infile, outfiles, pathname):
    for f in outfiles:
        os.unlink(f)
    outdir = os.path.dirname(pathname)
    outpattern = os.path.join(outdir, 'cluster_')
    minsl = globals_.config.get('Parameters', 'min_length_ratio')
    cmd = ('{usearch} -cluster_fast {infile} -id {id}'
           ' -clusters {outpattern} -sort length'
           ' -maxaccepts {max_accepts} -maxrejects {max_rejects} -minsl {minsl}')
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, id=globals_.config.get('Parameters', 'cluster_identity'),
                     outpattern=outpattern,
                     max_accepts=globals_.config.get('Parameters', 'max_accepts'),
                     max_rejects=globals_.config.get('Parameters', 'max_rejects'),
                     minsl=minsl)
    result = maybe_qsub(cmd, infile, outfiles, name="cluster")
    r = re.compile(r'^cluster_[0-9]+$')
    for f in list(f for f in os.listdir(outdir) if r.match(f)):
        oldfile = os.path.join(outdir, f)
        newfile = ''.join([oldfile, '.raw.fasta'])
        os.rename(oldfile, newfile)
    return result


@report_wrapper
def select_clusters(infile, outfile):
    records = list(SeqIO.parse(infile, 'fasta'))
    minsize = int(globals_.config.get('Parameters', 'min_cluster_size'))
    maxsize = int(globals_.config.get('Parameters', 'max_cluster_size'))
    if len(records) < minsize:
        touch(outfile)  # empty file; needed for pipeline sanity
        return
    if len(records) > maxsize:
        records = sample(records, maxsize)
    SeqIO.write(records, outfile, 'fasta')


# making wrappers like this is necessary because nested function
# definitions are not picklable.

@must_work(maybe=True)
@report_wrapper
def mafft_wrapper_maybe(infile, outfile):
    mafft(infile, outfile)


@must_work(seq_ids=True)
@report_wrapper
def mafft_wrapper_seq_ids(infile, outfile):
    mafft(infile, outfile)


@must_work(seq_ids=True)
@report_wrapper
def cat_wrapper_ids(infiles, outfile):
    cat_files(infiles, outfile)


@report_wrapper
def cat_wrapper(infiles, outfile):
    cat_files(infiles, outfile)


@must_work()
@report_wrapper
def translate_wrapper(infile, outfile):
    return translate_helper(infile, outfile, gapped=False, name='translate')


@must_work()
@report_wrapper
def gapped_translate_wrapper(infile, outfile):
    return translate_helper(infile, outfile, gapped=True, name='translate-gapped')


@must_work(maybe=True, illegal_chars='-')
@report_wrapper
def cluster_consensus(infile, outfile):
    n = re.search("cluster_([0-9]+)", infile).group(1)
    seq_id = next(SeqIO.parse(infile, 'fasta')).id
    label = re.split("_[0-9]+$", seq_id)[0]
    kwargs = {
        'python': globals_.config.get('Paths', 'python'),
        'script': os.path.join(globals_.script_dir, "DNAcons.py"),
        'infile': infile,
        'outfile': outfile,
        'ambifile': '{}.info'.format(outfile[:-len('.fasta')]),
        'id_str': "{label}_hqcs_{n}".format(label=label, n=n),
        }
    cmd = ("{python} {script} -o {outfile} --ambifile {ambifile}"
           " --id {id_str} {infile}".format(**kwargs))
    return maybe_qsub(cmd, infile, outfile, name="cluster-consensus")


def is_in_frame(seq):
    if len(seq) % 3 != 0:
        return False
    t = seq.translate()
    if len(t) != len(seq) / 3:
        return False
    if '*' in t:
        return False
    return True


def new_record_id(record, new_id):
    result = record[:]
    result.id = new_id
    result.name = ""
    result.description = ""
    return result


@must_work()
@report_wrapper
def inframe_hqcs(infile, outfile):
    """Filter to in-frame hqcs sequences with no stop codons"""
    records = SeqIO.parse(infile, 'fasta')
    result = list(new_record_id(r, '{}_inframe'.format(r.id))
                  for r in records if is_in_frame(r.seq))
    SeqIO.write(result, outfile, 'fasta')


def pause_for_editing_inframe_hqcs():
    infile = os.path.join(pipeline_dir, 'inframe_hqcs.fasta')
    outfile = os.path.join(pipeline_dir, 'inframe_hqcs.edited.fasta')
    if globals_.config.getboolean('Tasks', 'use_inframe_hqcs'):
        input('Paused for manual editing of {}'
              '\nPress Enter to continue.'.format(outfile))
        # check that result is not eempty
        n = sum(1 for r in SeqIO.parse(outfile, 'fasta'))
        if n == 0:
            raise Exception('{} is empty'.format(outfile))
    else:
        # just use global reference file
        dbfile = globals_.config.get('Parameters', 'reference_db')
        shutil.copyfile(dbfile, outfile)


@must_work()
@report_wrapper
def hqcs_db_search(infiles, outfile):
    infile, dbfile = infiles
    check_basename(dbfile, 'inframe_hqcs.edited.fasta')
    return usearch_reference_db(infile, outfile, dbfile=dbfile, name='hqcs-db-search')


@must_work(in_frame=True, illegal_chars='-')
@report_wrapper
def hqcs_shift_correction(infile, outfile):
    return shift_correction_helper(infile, outfile, keep=False, correct=True,
                                   name="hqcs-shift-correction")


@must_work(min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
@report_wrapper
def unique_hqcs(infile, outfile):
    cmd = ('{usearch} -derep_fulllength {infile} -fastaout {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    return maybe_qsub(cmd, infile, outfiles=outfile, name='unique_hqcs')


@must_work()
@report_wrapper
def make_individual_dbs(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    return maybe_qsub(cmd, infile, outfiles=outfile, name='make-individual-dbs')


def usearch_hqcs_ids(infile, outfile, dbfile, name=None):
    """run usearch_global against a hqcs database and print pairs of ids"""
    identity = globals_.config.get('Parameters', 'ccs_to_hqcs_identity')
    max_accepts = globals_.config.get('Parameters', 'max_accepts')
    max_rejects = globals_.config.get('Parameters', 'max_rejects')
    maxqt = globals_.config.get('Parameters', 'max_query_target_length_ratio')
    cmd = ("{usearch} -usearch_global {infile} -db {db} -id {id}"
           " -userout {outfile} -top_hit_only -userfields query+target -strand both"
           " -maxaccepts {max_accepts} -maxrejects {max_rejects}"
           " -maxqt {maxqt}")
    if name is None:
        name = 'usearch-hqcs-ids'
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, db=dbfile, id=identity, outfile=outfile,
                     max_accepts=max_accepts, max_rejects=max_rejects,
                     maxqt=maxqt)
    return maybe_qsub(cmd, infile, outfiles=outfile, name=name)


# FIXME: this is run with remote job limiter, but part of its task is run locally
@must_work()
@report_wrapper
def compute_copynumbers(infiles, outfile):
    hqcsfile, dbfile, ccsfile = infiles
    # make sure this suffix changes depending on what task comes before
    check_suffix(hqcsfile, '.uniques.fasta')
    check_suffix(dbfile, '.udb')
    check_suffix(ccsfile, '.lenfilter.fasta')
    pairfile = '{}.pairs'.format(remove_suffix(outfile, '.tsv'))
    usearch_hqcs_ids(ccsfile, pairfile, dbfile, name='compute-copynumber')
    with open(pairfile) as f:
        pairs = list(line.strip().split("\t") for line in f.readlines())
    hqcs_counts = defaultdict(lambda: 0)
    for pair in pairs:
        if len(pair) != 2:
            warnings.warn('CCS {} did not match any HQCS'.format(pair))
            continue
        ccs_id, hqcs_id = pair
        hqcs_counts[hqcs_id] += 1
    # deal with hqcs sequences with no copynumber by giving them 0
    ids = list(r.id for r in SeqIO.parse(hqcsfile, 'fasta'))
    for i in ids:
        if i not in hqcs_counts:
            hqcs_counts[i] = 0
    with open(outfile, 'w') as handle:
        for id_, count in hqcs_counts.items():
            handle.write('{}\t{}\n'.format(id_, count))


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
    stdout = os.path.join(globals_.qsub_dir, 'backtranslate.stdout')
    stderr = os.path.join(globals_.qsub_dir, 'backtranslate.stderr')
    return maybe_qsub(cmd, infiles, outfiles=outfile,
                      stdout=stdout, stderr=stderr,
                      name="backtranslate")


@must_work()
@report_wrapper
def make_hqcs_db(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    return maybe_qsub(cmd, infile, outfiles=outfile, name='make_hqcs_db')


@must_work()
@report_wrapper
def hqcs_ccs_pairs(infiles, outfile):
    # FIXME: this does basically the same thing as the copynumber task,
    # except it allows CCSs to map to HQCSs in different timepoints
    infile, dbfile = infiles
    check_suffix(infile, '.lenfilter.fasta')
    check_suffix(dbfile, '.degapped.udb')
    usearch_hqcs_ids(infile, outfile, dbfile, name='hqcs_ccs_ids')


@must_work()
@report_wrapper
def combine_pairs(infiles, outfiles, outdir):
    # FIXME: code duplication with compute_copynumbers
    for f in outfiles:
        os.unlink(f)
    pairfile, hqcsfile = infiles
    hqcs_suffix = '.hqcs-ccs-pairs.txt'
    ccsfile = "{}.fasta".format(pairfile[:-len(hqcs_suffix)])
    check_suffix(pairfile, hqcs_suffix)
    check_suffix(ccsfile, '.lenfilter.fasta')
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
    ccs_records = list(SeqIO.parse(ccsfile, "fasta"))
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
    infile, backtranslated = infiles
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
    hqcs, ccs, cn = infiles
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


def make_alignment_pipeline(name=None):
    if name is None:
        name = "alignment_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    filter_fastq_task = pipeline.transform(filter_fastq,
                                           input=None,
                                           filter=formatter(),
                                           output=os.path.join(pipeline_dir, '{basename[0]}.qfilter.fastq'))
    filter_fastq_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    trim_task = pipeline.transform(trim_heads_and_tails,
                                    input=filter_fastq_task,
                                    filter=suffix('.fastq'),
                                    output='.trimmed.fastq')
    trim_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    filter_runs_task = pipeline.transform(filter_runs,
                                          input=trim_task,
                                          filter=suffix('.fastq'),
                                          output='.noruns.fastq')
    filter_runs_task.jobs_limit(n_local_jobs, local_job_limiter)

    filter_contaminants_task = pipeline.transform(filter_contaminants,
                                                  input=filter_runs_task,
                                                  filter=suffix('.fasta'),
                                                  output='.good.fasta')
    filter_contaminants_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    filter_uncontaminated_task = pipeline.transform(filter_uncontaminated,
                                                    input=filter_contaminants_task,
                                                    filter=suffix('.fasta'),
                                                    output='.refpairs.fasta')
    filter_uncontaminated_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    if globals_.config.getboolean('Tasks', 'shift_correct_ccs'):
        ccs_shift_correction_task = pipeline.transform(ccs_shift_correction,
                                                       input=filter_uncontaminated_task,
                                                       filter=suffix('.fasta'),
                                                       output='.shift.fasta')
        ccs_shift_correction_task.jobs_limit(n_remote_jobs, remote_job_limiter)
        degap_input_task = ccs_shift_correction_task
    else:
        every_other_task = pipeline.transform(every_other,
                                              input=filter_uncontaminated_task,
                                              filter=suffix('.fasta'),
                                              output='.norefs.fasta')
        every_other_task.jobs_limit(n_local_jobs, local_job_limiter)
        degap_input_task = every_other_task

    degap_task = pipeline.transform(degap,
                                    input = degap_input_task,
                                    filter=suffix('.fasta'),
                                    output='degapped.fasta')
    degap_task.jobs_limit(n_local_jobs, local_job_limiter)

    filter_length_task = pipeline.transform(filter_length,
                                            input=degap_task,
                                            filter=suffix('.fasta'),
                                            output='.lenfilter.fasta')
    filter_length_task.jobs_limit(n_local_jobs, local_job_limiter)

    cluster_task = pipeline.subdivide(cluster,
                                      input=filter_length_task,
                                      filter=formatter('.*/(?P<LABEL>.+).qfilter'),
                                      output='{path[0]}/{LABEL[0]}.clusters/cluster_*.raw.fasta',
                                      extras=['{path[0]}/{LABEL[0]}.clusters/cluster_*.raw.fasta'])
    cluster_task.jobs_limit(n_remote_jobs, remote_job_limiter)
    cluster_task.mkdir(filter_length_task, formatter('.*/(?P<LABEL>.+).qfilter'), '{path[0]}/{LABEL[0]}.clusters')

    select_clusters_task = pipeline.transform(select_clusters,
                                              input=cluster_task,
                                              filter=suffix('.fasta'),
                                              output='.keep.fasta')
    select_clusters_task.jobs_limit(n_local_jobs, local_job_limiter)

    align_clusters_task = pipeline.transform(mafft_wrapper_maybe,
                                             name="align_clusters",
                                             input=select_clusters_task,
                                             filter=suffix('.fasta'),
                                             output='.aligned.fasta')
    align_clusters_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    cluster_consensus_task = pipeline.transform(cluster_consensus,
                                                input=align_clusters_task,
                                                filter=suffix('.fasta'),
                                                output='.hqcs.fasta')
    cluster_consensus_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    cat_clusters_task = pipeline.collate(cat_wrapper_ids,
                                         name="cat_clusters",
                                         input=cluster_consensus_task,
                                         filter=formatter(),
                                         output='{path[0]}.hqcs.fasta')
    cat_clusters_task.jobs_limit(n_local_jobs, local_job_limiter)

    inframe_hqcs_task = pipeline.transform(inframe_hqcs,
                                           input=cat_clusters_task,
                                           filter=suffix('.fasta'),
                                           output='.inframe.fasta')
    inframe_hqcs_task.jobs_limit(n_local_jobs, local_job_limiter)

    cat_inframe_hqcs_task = pipeline.merge(cat_wrapper_ids,
                                           input=inframe_hqcs_task,
                                           output=os.path.join(pipeline_dir, 'inframe_hqcs.fasta'))
    cat_inframe_hqcs_task.jobs_limit(n_local_jobs, local_job_limiter)

    copy_inframe_hqcs_task = pipeline.transform(copy_file,
                                                input=cat_inframe_hqcs_task,
                                                filter=suffix('.fasta'),
                                                output='.edited.fasta')
    copy_inframe_hqcs_task.jobs_limit(n_local_jobs, local_job_limiter)
    copy_inframe_hqcs_task.posttask(pause_for_editing_inframe_hqcs)

    hqcs_db_search_task = pipeline.transform(hqcs_db_search,
                                             input=cat_clusters_task,
                                             add_inputs=add_inputs(copy_inframe_hqcs_task),
                                             filter=suffix('.fasta'),
                                             output='.refpairs.fasta')
    hqcs_db_search_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    hqcs_shift_correction_task = pipeline.transform(hqcs_shift_correction,
                                                    input=hqcs_db_search_task,
                                                    filter=suffix('.fasta'),
                                                    output='.shifted.fasta')
    hqcs_shift_correction_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    unique_hqcs_task = pipeline.transform(unique_hqcs,
                                          input=hqcs_shift_correction_task,
                                          filter=suffix('.fasta'),
                                          output='.uniques.fasta')
    unique_hqcs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    make_individual_dbs_task = pipeline.transform(make_individual_dbs,
                                                  input=unique_hqcs_task,
                                                  filter=suffix('.fasta'),
                                                  output='.udb')
    make_individual_dbs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    compute_copynumbers_task = pipeline.collate(compute_copynumbers,
                                                input=[unique_hqcs_task, make_individual_dbs_task, filter_length_task],
                                                filter=formatter(r'.*/(?P<LABEL>.+).(qfilter|clusters)'),
                                                output=os.path.join(pipeline_dir, '{LABEL[0]}.copynumbers.tsv'))
    compute_copynumbers_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    merge_copynumbers_task = pipeline.merge(cat_wrapper,
                                            name='cat_copynumbers',
                                            input=compute_copynumbers_task,
                                            output=os.path.join(pipeline_dir, 'copynumbers.tsv'))
    merge_copynumbers_task.jobs_limit(n_local_jobs, local_job_limiter)

    cat_all_hqcs_task = pipeline.merge(cat_all_hqcs,
                                       input=unique_hqcs_task,
                                       output=os.path.join(pipeline_dir, "hqcs.fasta"))
    cat_all_hqcs_task.jobs_limit(n_local_jobs, local_job_limiter)

    translate_hqcs_task = pipeline.transform(translate_wrapper,
                                             name='translate_hqcs',
                                             input=cat_all_hqcs_task,
                                             filter=suffix('.fasta'),
                                             output='.translated.fasta')
    translate_hqcs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    align_hqcs_protein_task = pipeline.transform(mafft_wrapper_seq_ids,
                                                 name='align_hqcs_protein',
                                                 input=translate_hqcs_task,
                                                 filter=suffix('.fasta'),
                                                 output='.aligned.fasta')
    align_hqcs_protein_task.jobs_limit(n_local_jobs, local_job_limiter)

    copy_protein_alignment_task = pipeline.transform(copy_file,
                                                     name='copy_protein_alignment',
                                                     input=align_hqcs_protein_task,
                                                     filter=suffix('.fasta'),
                                                     output='.edited.fasta')
    copy_protein_alignment_task.jobs_limit(n_local_jobs, local_job_limiter)
    copy_protein_alignment_task.posttask(pause_for_editing_alignment)

    backtranslate_alignment_task = pipeline.merge(backtranslate_alignment,
                                                  input=[cat_all_hqcs_task,
                                                         copy_protein_alignment_task],
                                                  output=os.path.join(pipeline_dir, 'hqcs.translated.aligned.edited.backtranslated.fasta'))
    backtranslate_alignment_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    pipeline.set_head_tasks([filter_fastq_task])
    pipeline.set_tail_tasks([backtranslate_alignment_task, merge_copynumbers_task])

    if globals_.config.getboolean('Tasks', 'align_ccs'):
        degap_backtranslated_alignment_task = pipeline.transform(degap,
                                                                 input=backtranslate_alignment_task,
                                                                 filter=suffix('.fasta'),
                                                                 output='.degapped.fasta')
        degap_backtranslated_alignment_task.jobs_limit(n_local_jobs, local_job_limiter)

        make_hqcs_db_task = pipeline.transform(make_hqcs_db,
                                               input=degap_backtranslated_alignment_task,
                                               filter=suffix('.fasta'),
                                               output='.udb')
        make_hqcs_db_task.jobs_limit(n_remote_jobs, remote_job_limiter)

        hqcs_ccs_pairs_task = pipeline.transform(hqcs_ccs_pairs,
                                                 input=filter_length_task,
                                                 filter=suffix('.fasta'),
                                                 add_inputs=add_inputs(make_hqcs_db),
                                                 output='.hqcs-ccs-pairs.txt')
        hqcs_ccs_pairs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

        combine_pairs_task = pipeline.subdivide(combine_pairs,
                                                input=[hqcs_ccs_pairs_task],
                                                filter=formatter('.*/(?P<LABEL>.+).qfilter'),
                                                add_inputs=add_inputs(degap_backtranslated_alignment),
                                                output='{path[0]}/{LABEL[0]}.ccs-alignments/combined.*.unaligned.fasta',
                                                extras=['{path[0]}/{LABEL[0]}.ccs-alignments'])
        combine_pairs_task.jobs_limit(n_local_jobs, local_job_limiter)
        combine_pairs_task.mkdir(hqcs_ccs_pairs_task,
                                 formatter('.*/(?P<LABEL>.+).qfilter'),
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

        insert_gaps_task = pipeline.transform(insert_gaps_wrapper,
                                              input=convert_bam_to_fasta_task,
                                              filter=suffix('.fasta'),
                                              add_inputs=add_inputs(backtranslate_alignment),
                                              output='.gapped.fasta')
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

        diagnosis_output = list(os.path.join(pipeline_dir,
                                             'diagnosis',
                                             "freq_agreement_no_x_{}.png".format(t.label))
                                for t in globals_.timepoints)
        diagnose_alignment_task = pipeline.merge(diagnose_alignment,
                                                 name='diagnose_alignment',
                                                 input=[copy_protein_alignment_task,
                                                        translate_ccs_task,
                                                        merge_copynumbers_task],
                                                 output=diagnosis_output)
        diagnose_alignment_task.jobs_limit(n_remote_jobs, remote_job_limiter)
        diagnose_alignment_task.mkdir(os.path.join(pipeline_dir, 'diagnosis'))

    for task in pipeline.head_tasks:
        task.mkdir(pipeline_dir)

    return pipeline
