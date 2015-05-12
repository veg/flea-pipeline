import os
import re
from random import sample
import tempfile
from functools import partial
from collections import defaultdict
import json

from ruffus import Pipeline, suffix, formatter, add_inputs

from Bio import SeqIO

from util import maybe_qsub, cat_files, touch, strlist, traverse
from util import new_record_seq_str, insert_gaps, update_record_seq
from util import must_work, must_produce
from util import local_job_limiter, remote_job_limiter
from util import check_suffix, check_basename, n_jobs

import pipeline_globals as globals_

from translate import translate
from backtranslate import backtranslate
from DNAcons import consfile
from correct_shifts import correct_shifts_fasta, write_correction_result
from perfectORFs import perfect_file


def mafft(infile, outfile):
    stderr = '{}.stderr'.format(outfile)
    cmd = 'mafft-fftns --ep 0.5 --quiet --preservecase {} > {} 2>{}'.format(infile, outfile, stderr)
    maybe_qsub(cmd, outfiles=outfile, stdout='/dev/null', stderr='/dev/null')


@must_work()
def filter_fastq(infile, outfile):
    qmax = globals_.config.get('Parameters', 'qmax')
    min_len = globals_.config.get('Parameters', 'min_sequence_length')
    max_err_rate = globals_.config.get('Parameters', 'max_err_rate')
    cmd = ('{usearch} -fastq_filter {infile} -fastq_maxee_rate {max_err_rate}'
         ' -threads 1 -fastq_qmax {qmax} -fastq_minlen {min_len} -fastaout {outfile}'
         ' -relabel "{seq_id}_"'.format(
            usearch=globals_.config.get('Paths', 'usearch'),
            infile=infile, outfile=outfile, qmax=qmax, min_len=min_len,
            max_err_rate=max_err_rate, seq_id=globals_.timepoint_ids[infile]))
    maybe_qsub(cmd, outfiles=outfile, name='filter-fastq')


def filter_contaminants(infile, outfiles):
    uncontam, contam = outfiles
    cmd = ('{usearch} -usearch_global {infile} -db {db}'
           ' -id {id} -notmatched {uncontam} -matched {contam}'
           ' -strand both'.format(usearch=globals_.config.get('Paths', 'usearch'),
                                  infile=infile, db=globals_.config.get('Parameters', 'contaminants_db'),
                                  id=globals_.config.get('Parameters', 'contaminant_identity'),
                                  uncontam=uncontam, contam=contam))
    maybe_qsub(cmd, outfiles=outfiles, name='filter-contaminants')


def usearch_global_pairs(infile, outfile, dbfile, identity, nums_only=False, name=None):
    max_accepts = globals_.config.get('Parameters', 'max_accepts')
    max_rejects = globals_.config.get('Parameters', 'max_rejects')
    if nums_only:
        cmd = ("{usearch} -usearch_global {infile} -db {db} -id {id}"
               " -userout {outfile} -top_hit_only -userfields query+target -strand both"
               " -maxaccepts {max_accepts} -maxrejects {max_rejects}")
    else:
        cmd = ('{usearch} -usearch_global {infile} -db {db} -id {id}'
               ' -fastapairs {outfile} -top_hit_only -strand both'
               ' -maxaccepts {max_accepts} -maxrejects {max_rejects}')
    if name is None:
        name = 'usearch-global-pairs'
    cmd = cmd.format(usearch=globals_.config.get('Paths', 'usearch'),
                     infile=infile, db=dbfile, id=identity, outfile=outfile,
                     max_accepts=max_accepts, max_rejects=max_rejects)
    maybe_qsub(cmd, outfiles=outfile, name=name)


def usearch_global_get_pairs(infile, dbfile, identity, name=None):
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, 'pairs.txt')
        usearch_global_pairs(infile, outfile, dbfile, identity, nums_only=True, name=name)
        with open(outfile) as f:
            return list(line.strip().split("\t") for line in f.readlines())


@must_work()
def filter_uncontaminated(infiles, outfile):
    uncontam, _ = infiles
    usearch_global_pairs(uncontam, outfile, globals_.config.get('Parameters', 'reference_db'),
                         globals_.config.get('Parameters', 'reference_identity'),
                         name="filter-uncontaminated")


@must_work(seq_ratio=(2, 1))
def shift_correction(infile, outfile):
    correct_shifts_fasta(infile, outfile, keep=True)


@must_work(seq_ids=True)
def sort_by_length(infile, outfile):
    cmd = ('{usearch} -sortbylength {infile} -fastaout {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    maybe_qsub(cmd, outfiles=outfile, name="sort-by-length")


@must_produce(n=int(globals_.config.get('Parameters', 'min_n_clusters')))
def cluster(infile, outfiles, pathname):
    for f in outfiles:
        os.unlink(f)
    outdir = '{}.clusters'.format(infile[:-len('.fasta')])
    outpattern = os.path.join(outdir, 'cluster_')
    cmd = ('{usearch} -cluster_fast {infile} -id {id}'
           ' -clusters {outpattern} -maxaccepts {max_accepts}'
           ' -maxrejects {max_rejects}'.format(
            usearch=globals_.config.get('Paths', 'usearch'),
            infile=infile, id=globals_.config.get('Parameters', 'cluster_identity'),
            outpattern=outpattern,
            max_accepts=globals_.config.get('Parameters', 'max_accepts'),
            max_rejects=globals_.config.get('Parameters', 'max_rejects')))
    maybe_qsub(cmd, name="cluster")
    r = re.compile(r'^cluster_[0-9]+$')
    for f in list(f for f in os.listdir(outdir) if r.match(f)):
        oldfile = os.path.join(outdir, f)
        newfile = ''.join([oldfile, '.raw.fasta'])
        os.rename(oldfile, newfile)


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
def mafft_wrapper_maybe(infile, outfile):
    mafft(infile, outfile)


@must_work(seq_ids=True)
def mafft_wrapper_seq_ids(infile, outfile):
    mafft(infile, outfile)


@must_work(seq_ids=True)
def cat_wrapper(infiles, outfile):
    cat_files(infiles, outfile)


@must_work()
def translate_wrapper(infile, outfile):
    translate(infile, outfile)


@must_work(maybe=True, illegal_chars='-')
def cluster_consensus(infile, outfile):
    ambifile = '{}.info'.format(outfile[:-len('.fasta')])
    consfile(infile, outfile, ambifile, ungap=True)


@must_work()
def consensus_db_search(infile, outfile):
    usearch_global_pairs(infile, outfile, globals_.config.get('Parameters', 'reference_db'),
                         globals_.config.get('Parameters', 'reference_identity'),
                         name='consensus-db-search')


@must_work()
def consensus_shift_correction(infile, outfile):
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile, keep=False)
    sumfile = '{}.summary'.format(outfile)
    write_correction_result(n_seqs, n_fixed, sumfile)


@must_work(seq_ids=True)
def sort_consensus(infile, outfile):
    cmd = ('{usearch} -sortbylength {infile} -fastaout {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    maybe_qsub(cmd, outfiles=outfile, name='sort-consensus')


@must_work(min_seqs=int(globals_.config.get('Parameters', 'min_n_clusters')))
def unique_consensus(infile, outfile):
    cmd = ('{usearch} -cluster_fast {infile} -id 1 -centroids {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    maybe_qsub(cmd, outfiles=outfile, name='unique_consensus')


@must_work()
def perfect_orfs(infile, outfile):
    perfect_file(infile, outfile, min_len=int(globals_.config.get('Parameters', 'min_orf_length')),
                 table=1, verbose=False)


@must_work()
def make_individual_dbs(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    maybe_qsub(cmd,  outfiles=outfile, name='make-individual-dbs')


# FIXME: this is run with remote job limiter, but part of its task is run locally
@must_work()
def compute_copynumbers(infiles, outfile, basename):
    rawfile, perfectfile, dbfile = infiles
    check_suffix(perfectfile, '.perfect.fasta')
    check_suffix(dbfile, '.udb')
    identity = globals_.config.get('Parameters', 'raw_to_consensus_identity')
    pairfile = '{}.copynumber.pairs'.format(basename)
    usearch_global_pairs(rawfile, pairfile, dbfile, identity,
                         nums_only=True, name='compute-copynumber')
    with open(pairfile) as f:
            pairs = list(line.strip().split("\t") for line in f.readlines())
    # FIXME: what if no sequence maps to a consensus?
    consensus_counts = defaultdict(lambda: 0)
    for raw_id, ref_id in pairs:
        consensus_counts[ref_id] += 1
    with open(outfile, 'w') as handle:
        json.dump(consensus_counts, handle, separators=(",\n", ":"))


@must_work()
def merge_copynumbers(infiles, outfile):
    result = {}
    for infile in infiles:
        with open(infile) as handle:
            result.update(json.load(handle))
    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


@must_work(seq_ids=True)
def cat_all_perfect(infiles, outfile):
    if len(infiles) != len(globals_.timepoints):
        raise Exception('Number of input files does not match number'
                        ' of timepoints')
    cat_files(infiles, outfile)


def pause(filename):
    if globals_.config.getboolean('Tasks', 'pause_after_codon_alignment'):
        input('Paused for manual editing of {}'
              '\nPress Enter to continue.'.format(filename))


@must_work()
def backtranslate_alignment(infiles, outfile):
    perfect, aligned = infiles
    check_basename(perfect, 'all_perfect_orfs.fasta')
    backtranslate(aligned, perfect, outfile)


@must_work()
def degap_backtranslated_alignment(infile, outfile):
    # we need this in case sequences were removed during hand-editing
    # of the output of `codon_align_perfect()`.
    records = (SeqIO.parse(infile, 'fasta'))
    processed = (update_record_seq(r, r.seq.ungap('-')) for r in records)
    SeqIO.write(processed, outfile, 'fasta')


@must_work()
def make_full_db(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    maybe_qsub(cmd, outfiles=outfile, name='make_full_db')



@must_work()
def full_timestep_pairs(infiles, outfile):
    infile, dbfile = infiles
    identity = globals_.config.get('Parameters', 'raw_to_consensus_identity')
    usearch_global_pairs(infile, outfile, dbfile, identity, nums_only=True,
                         name='full-timestep-pairs')


@must_work()
def combine_pairs(infiles, outfiles, basename):
    for f in outfiles:
        os.unlink(f)
    infile, perfectfile = infiles
    seqsfile = '{}.fasta'.format(basename)
    with open(infile) as handle:
        pairs = list(line.strip().split("\t") for line in handle.readlines())
    match_dict = defaultdict(list)
    for seq, ref in pairs:
        match_dict[ref].append(seq)
    references_records = list(SeqIO.parse(perfectfile, "fasta"))
    seq_records = list(SeqIO.parse(seqsfile, "fasta"))
    references_dict = {r.id : r for r in references_records}
    seq_dict = {r.id : r for r in seq_records}

    for ref_id, seq_ids in match_dict.items():
        outfile = '{}.alignments/combined.{}.unaligned.fasta'.format(basename, ref_id)
        records = [references_dict[ref_id]]
        records.extend(list(seq_dict[i] for i in seq_ids))
        SeqIO.write(records, outfile, "fasta")


@must_work()
def codon_align(infile, outfile):
    cmd = "{} -R {} {}".format(globals_.config.get('Paths', 'bealign'), infile, outfile)
    stdout = os.path.join(globals_.qsub_dir, '{}.stdout'.format(outfile))
    stderr = os.path.join(globals_.qsub_dir, '{}.stderr'.format(outfile))
    maybe_qsub(cmd, outfiles=outfile, stdout=stdout, stderr=stderr)


@must_work()
def convert_bam_to_fasta(infile, outfile):
    cmd = "{} {} {}".format(globals_.config.get('Paths', 'bam2msa'), infile, outfile)
    stdout = os.path.join(globals_.qsub_dir, '{}.stdout'.format(outfile))
    stderr = os.path.join(globals_.qsub_dir, '{}.stderr'.format(outfile))
    maybe_qsub(cmd, outfiles=outfile, stdout=stdout, stderr=stderr)


@must_work()
def insert_gaps_wrapper(infiles, outfile):
    infile, backtranslated = infiles
    ref, *seqs = list(SeqIO.parse(infile, 'fasta'))
    ref_gapped = list(r for r in SeqIO.parse(backtranslated, 'fasta')
                      if r.id == ref.id)[0]
    seqs_gapped = (new_record_seq_str(r, insert_gaps(str(ref_gapped.seq),
                                                     str(r.seq),
                                                     '-', '-', skip=1))
                   for r in seqs)
    SeqIO.write(seqs_gapped, outfile, "fasta")


def make_alignment_pipeline(name=None):
    if name is None:
        name = "alignment_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    filter_fastq_task = pipeline.transform(filter_fastq,
                                           input=None,
                                           filter=suffix(".fastq"),
                                           output='.qfilter.fasta')
    filter_fastq_task.jobs_limit(n_local_jobs, local_job_limiter)
    pipeline.set_head_tasks([filter_fastq_task])

    filter_contaminants_task = pipeline.transform(filter_contaminants,
                                                  input=filter_fastq_task,
                                                  filter=suffix('.fasta'),
                                                  output=['.uncontam.fasta', '.contam.fasta'],)
    filter_contaminants_task.jobs_limit(n_remote_jobs, remote_job_limiter)


    filter_uncontaminated_task = pipeline.transform(filter_uncontaminated,
                                                    input=filter_contaminants_task,
                                                    filter=suffix('.uncontam.fasta'),
                                                    output='.uncontam.rfilter.fasta')
    filter_uncontaminated_task.jobs_limit(n_remote_jobs, remote_job_limiter)


    shift_correction_task = pipeline.transform(shift_correction,
                                               input=filter_uncontaminated_task,
                                               filter=suffix('.fasta'),
                                               output='.shifted.fasta')
    shift_correction_task.jobs_limit(n_local_jobs, local_job_limiter)


    sort_by_length_task = pipeline.transform(sort_by_length,
                                             input=shift_correction_task,
                                             filter=suffix('.fasta'),
                                             output='.sorted.fasta')
    sort_by_length_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    cluster_task = pipeline.subdivide(cluster,
                                      input=sort_by_length_task,
                                      filter=formatter(),
                                      output='{path[0]}/{basename[0]}.clusters/cluster_*.raw.fasta',
                                      extras=['{path[0]}/{basename[0]}.clusters/cluster_*.raw.fasta'])
    cluster_task.jobs_limit(n_remote_jobs, remote_job_limiter)
    cluster_task.mkdir(sort_by_length_task, suffix('.fasta'), '.clusters')


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
                                                output='.cons.fasta')
    cluster_consensus_task.jobs_limit(n_local_jobs, local_job_limiter)


    cat_clusters_task = pipeline.collate(cat_wrapper,
                                         name="cat_clusters",
                                         input=cluster_consensus_task,
                                         filter=formatter(),
                                         output='{path[0]}.cons.fasta')
    cat_clusters_task.jobs_limit(n_local_jobs, local_job_limiter)


    consensus_db_search_task = pipeline.transform(consensus_db_search,
                                                  input=cat_clusters_task,
                                                  filter=suffix('.fasta'),
                                                  output='.pairs.fasta')
    consensus_db_search_task.jobs_limit(n_remote_jobs, remote_job_limiter)


    consensus_shift_correction_task = pipeline.transform(consensus_shift_correction,
                                                         input=consensus_db_search_task,
                                                         filter=suffix('.fasta'),
                                                         output='.shifted.fasta')
    consensus_shift_correction_task.jobs_limit(n_local_jobs, local_job_limiter)


    sort_consensus_task = pipeline.transform(sort_consensus,
                                             input=consensus_shift_correction_task,
                                             filter=suffix('.fasta'),
                                             output='.sorted.fasta')
    sort_consensus_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    unique_consensus_task = pipeline.transform(unique_consensus,
                                               input=sort_consensus_task,
                                               filter=suffix('.fasta'),
                                               output='.uniques.fasta')
    unique_consensus_task.jobs_limit(n_remote_jobs, remote_job_limiter)


    perfect_orfs_task = pipeline.transform(perfect_orfs,
                                           input=unique_consensus_task,
                                           filter=suffix('.fasta'),
                                           output='.perfect.fasta')
    perfect_orfs_task.jobs_limit(n_local_jobs, local_job_limiter)


    make_individual_dbs_task = pipeline.transform(make_individual_dbs,
                                                  input=perfect_orfs_task,
                                                  filter=suffix('.fasta'),
                                                  output='.udb')
    make_individual_dbs_task.jobs_limit(n_remote_jobs, remote_job_limiter)


    compute_copynumbers_task = pipeline.collate(compute_copynumbers,
                                               input=[perfect_orfs_task, make_individual_dbs_task, shift_correction_task],
                                               filter=formatter(r'(?P<NAME>.+).qfilter'),
                                               output='{NAME[0]}.copynumbers.json',
                                               extras=['{NAME[0]}'])
    compute_copynumbers_task.jobs_limit(n_remote_jobs, remote_job_limiter)

    merge_copynumbers_task = pipeline.merge(merge_copynumbers,
                                            input=compute_copynumbers_task,
                                            output=os.path.join(globals_.data_dir, 'copynumbers.json'))
    merge_copynumbers_task.jobs_limit(n_local_jobs, local_job_limiter)

    cat_all_perfect_task = pipeline.merge(cat_all_perfect,
                                          input=perfect_orfs_task,
                                          output=os.path.join(globals_.data_dir, "all_perfect_orfs.fasta"))
    cat_all_perfect_task.jobs_limit(n_local_jobs, local_job_limiter)

    translate_perfect_task = pipeline.transform(translate_wrapper,
                                                name='translate_perfect',
                                                input=cat_all_perfect_task,
                                                filter=suffix('.fasta'),
                                                output='.translated.fasta')
    translate_perfect_task.jobs_limit(n_local_jobs, local_job_limiter)

    codon_align_perfect_task = pipeline.transform(mafft_wrapper_seq_ids,
                                                  name='codon_align_perfect',
                                                  input=translate_perfect_task,
                                                  filter=suffix('.fasta'),
                                                  output='.aligned.fasta')
    codon_align_perfect_task.jobs_limit(n_local_jobs, local_job_limiter)
    codon_align_perfect_task.posttask(partial(pause, 'protein alignment'))

    backtranslate_alignment_task = pipeline.merge(backtranslate_alignment,
                                              input=[cat_all_perfect_task,
                                                     codon_align_perfect_task],
                                              output=os.path.join(globals_.data_dir, 'all_backtranslated.fas'))
    backtranslate_alignment_task.jobs_limit(n_local_jobs, local_job_limiter)

    pipeline.set_tail_tasks([backtranslate_alignment_task])

    if globals_.config.getboolean('Tasks', 'align_full'):
        degap_backtranslated_alignment_task = pipeline.transform(degap_backtranslated_alignment,
                                                                 input=backtranslate_alignment_task,
                                                                 filter=formatter(),
                                                                 output='{path[0]}/all_perfect_orfs_degapped.fasta')
        degap_backtranslated_alignment_task.jobs_limit(n_local_jobs, local_job_limiter)

        make_full_db_task = pipeline.transform(make_full_db,
                                               input=degap_backtranslated_alignment_task,
                                               filter=suffix('.fasta'),
                                               output='.udb')
        make_full_db_task.jobs_limit(n_remote_jobs, remote_job_limiter)


        full_timestep_pairs_task = pipeline.transform(full_timestep_pairs,
                                                      input=shift_correction_task,
                                                      filter=suffix('.fasta'),
                                                      add_inputs=add_inputs(make_full_db),
                                                      output='.pairs.txt')
        full_timestep_pairs_task.jobs_limit(n_remote_jobs, remote_job_limiter)

        combine_pairs_task = pipeline.subdivide(combine_pairs,
                                                input=full_timestep_pairs_task,
                                                filter=formatter('.*/(?P<NAME>.+).pairs.txt'),
                                                add_inputs=add_inputs(degap_backtranslated_alignment),
                                                output='{path[0]}/{NAME[0]}.alignments/combined.*.unaligned.fasta',
                                                extras=['{path[0]}/{NAME[0]}'])
        combine_pairs_task.jobs_limit(n_local_jobs, local_job_limiter)
        combine_pairs_task.mkdir(full_timestep_pairs_task, suffix('.pairs.txt'), '.alignments')

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

        merge_all_timepoints_task = pipeline.merge(cat_wrapper,
                                                   name='merge_all_timepoints',
                                                   input=insert_gaps_task,
                                                   output=os.path.join(globals_.data_dir, 'all_timepoints.aligned.fasta'))
        merge_all_timepoints_task.jobs_limit(n_local_jobs, local_job_limiter)

    return pipeline
