import uuid
import shutil
import re
from random import sample


from fnmatch import fnmatch
import tempfile
from functools import partial

from Bio import SeqIO

from util import call, qsub, cat_files, touch, strlist, traverse
from util import new_record_seq_str, insert_gaps, update_record_seq
import util

import pipeline_globals as globals_

from translate import translate
from backtranslate import backtranslate
from DNAcons import consfile
from correct_shifts import correct_shifts_fasta, write_correction_result
from perfectORFs import perfect_file


def mafft(infile, outfile):
    stderr = '{}.stderr'.format(outfile)
    cmd = 'mafft-fftns --ep 0.5 --quiet --preservecase {} > {} 2>{}'.format(infile, outfile, stderr)
    maybe_qsub(cmd, globals_.config, outfiles=outfile, stdout='/dev/null', stderr='/dev/null')


@must_work()
def filter_fastq(infiles, outfile):
    infile, _ = infiles
    qmax = globals_.config.get('Parameters', 'qmax')
    min_len = globals_.config.get('Parameters', 'min_sequence_length')
    max_err_rate = globals_.config.get('Parameters', 'max_err_rate')
    cmd = ('{usearch} -fastq_filter {infile} -fastq_maxee_rate {max_err_rate}'
         ' -threads 1 -fastq_qmax {qmax} -fastq_minlen {min_len} -fastaout {outfile}'
         ' -relabel "{seq_id}_"'.format(
            usearch=globals_.config.get('Paths', 'usearch'),
            infile=infile, outfile=outfile, qmax=qmax, min_len=min_len,
            max_err_rate=max_err_rate, seq_id=timepoint_ids[infile]))
    maybe_qsub(cmd, globals_.config, outfiles=outfile, name='filter-fastq')


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


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(filter_contaminants, suffix('.uncontam.fasta'),
           '.uncontam.rfilter.fasta')
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
def cluster(infile, outfiles):
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


@must_work(maybe=True)
def align_clusters(infile, outfile):
    mafft(infile, outfile)


@must_work(maybe=True, illegal_chars='-')
def cluster_consensus(infile, outfile):
    ambifile = '{}.info'.format(outfile[:-len('.fasta')])
    consfile(infile, outfile, ambifile, ungap=True)


@must_work(seq_ids=True)
def cat_clusters(infiles, outfile):
    cat_files(infiles, outfile)


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


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(sort_consensus, suffix('.fasta'), '.uniques.fasta')
@must_work()
def unique_consensus(infile, outfile):
    cmd = ('{usearch} -cluster_fast {infile} -id 1 -centroids {outfile}'.format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    maybe_qsub(cmd, outfiles=outfile, name='unique_consensus')


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(unique_consensus, suffix('.fasta'), '.perfect.fasta')
@must_work()
def perfect_orfs(infile, outfile):
    perfect_file(infile, outfile, min_len=int(globals_.config.get('Parameters', 'min_orf_length')),
                 table=1, verbose=False)


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(perfect_orfs, suffix('.fasta'), '.udb')
@must_work()
def make_individual_dbs(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    maybe_qsub(cmd,  outfiles=outfile, name='make-individual-dbs')


@jobs_limit(n_remote_jobs, remote_job_limiter)
@collate([perfect_orfs, make_individual_dbs, shift_correction],
         formatter(r'(?P<NAME>.+).qfilter'),
         '{NAME[0]}.copynumber.fasta',
         '{NAME[0]}')
@must_work(seq_ratio=(1, 1), pattern='*perfect*fasta')
def add_copynumber(infiles, outfile, basename):
    rawfile, perfectfile, dbfile = infiles
    check_suffix(perfectfile, '.perfect.fasta')
    check_suffix(dbfile, '.udb')
    identity = globals_.config.get('Parameters', 'raw_to_consensus_identity')
    pairfile = '{}.copynumber.pairs'.format(basename)
    usearch_global_pairs(rawfile, pairfile, dbfile, identity,
                         nums_only=True, name='add-copynumber')
    with open(pairfile) as f:
            pairs = list(line.strip().split("\t") for line in f.readlines())
    consensus_counts = defaultdict(lambda: 1)
    for raw_id, ref_id in pairs:
        consensus_counts[str(ref_id).upper()] += 1
    def new_record(r):
        r = r[:]
        _id = str(r.id).upper()
        r.id = "{}_{}".format(_id, consensus_counts[_id])
        r.name = ''
        r.description = ''
        return r
    SeqIO.write((new_record(r) for r in SeqIO.parse(perfectfile, 'fasta')),
                outfile, 'fasta')


@jobs_limit(n_local_jobs, local_job_limiter)
@merge(add_copynumber, "all_perfect_orfs.fasta")
@must_work(seq_ids=True)
def cat_all_perfect(infiles, outfile):
    if len(infiles) != len(timepoints):
        raise Exception('Number of input files does not match number'
                        ' of timepoints')
    cat_files(infiles, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(cat_all_perfect, suffix('.fasta'), '.translated.fasta')
@must_work(seq_ids=True)
def translate_perfect(infile, outfile):
    translate(infile, outfile)


def pause(filename):
    if globals_.config.getboolean('Tasks', 'pause_after_codon_alignment'):
        input('Paused for manual editing of {}'
              '\nPress Enter to continue.'.format(filename))


@posttask(partial(pause, 'hyphy_data/input/merged.prot'))
@mkdir(hyphy_input_dir)
@mkdir(hyphy_results_dir)
@jobs_limit(n_local_jobs, local_job_limiter)
@transform(translate_perfect, formatter(), hyphy_input('merged.prot'))
@must_work(seq_ids=True)
def codon_align_perfect(infile, outfile):
    mafft(infile, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@merge([cat_all_perfect, codon_align_perfect], 'all_backtranslated.fas')
@must_work()
def backtranslate_alignment(infiles, outfile):
    perfect, aligned = infiles
    check_basename(perfect, 'all_perfect_orfs.fasta')
    backtranslate(aligned, perfect, outfile)


@jobs_limit(n_local_jobs, local_job_limiter)
@transform(backtranslate_alignment, formatter(), 'all_perfect_orfs_degapped.fasta')
@must_work()
def degap_backtranslated_alignment(infile, outfile):
    # we need this in case sequences were removed during hand-editing
    # of the output of `codon_align_perfect()`.
    records = (SeqIO.parse(infile, 'fasta'))
    processed = (update_record_seq(r, r.seq.ungap('-')) for r in records)
    SeqIO.write(processed, outfile, 'fasta')


@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(degap_backtranslated_alignment, suffix('.fasta'), '.udb')
@must_work()
def make_full_db(infile, outfile):
    cmd = ("{usearch} -makeudb_usearch {infile} -output {outfile}".format(
            usearch=globals_.config.get('Paths', 'usearch'), infile=infile, outfile=outfile))
    maybe_qsub(cmd, outfiles=outfile, name='make_full_db')



@active_if(globals_.config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(shift_correction,
           suffix('.fasta'),
           add_inputs([make_full_db]),
           '.pairs.txt')
@must_work()
def full_timestep_pairs(infiles, outfile):
    infile, (dbfile,) = infiles
    identity = globals_.config.get('Parameters', 'raw_to_consensus_identity')
    usearch_global_pairs(infile, outfile, dbfile, identity, nums_only=True,
                         name='full-timestep-pairs')


@active_if(globals_.config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_local_jobs, local_job_limiter)
@mkdir(full_timestep_pairs, suffix('.pairs.txt'), '.alignments')
@subdivide(full_timestep_pairs,
           formatter('.*/(?P<NAME>.+).pairs.txt'),
           add_inputs([degap_backtranslated_alignment]),
           '{NAME[0]}.alignments/combined.*.unaligned.fasta',
           '{NAME[0]}')
@must_work()
def combine_pairs(infiles, outfiles, basename):
    for f in outfiles:
        os.unlink(f)
    infile, (perfectfile,) = infiles
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


@active_if(globals_.config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(combine_pairs,
           suffix('.unaligned.fasta'),
           '.aligned.bam')
@must_work()
def codon_align(infile, outfile):
    cmd = "{} -R {} {}".format(globals_.config.get('Paths', 'bealign'), infile, outfile)
    stdout = '{}.stdout'.format(outfile)
    stderr = '{}.stderr'.format(outfile)
    maybe_qsub(cmd, outfiles=outfile, stdout=stdout, stderr=stderr)


@active_if(globals_.config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_remote_jobs, remote_job_limiter)
@transform(codon_align, suffix('.bam'), '.fasta')
@must_work()
def convert_bam_to_fasta(infile, outfile):
    cmd = "{} {} {}".format(globals_.config.get('Paths', 'bam2msa'), infile, outfile)
    stdout = '{}.stdout'.format(outfile)
    stderr = '{}.stderr'.format(outfile)
    maybe_qsub(cmd, outfiles=outfile, stdout=stdout, stderr=stderr)


@active_if(globals_.config.getboolean('Tasks', 'align_full'))
@jobs_limit(n_local_jobs, local_job_limiter)
@transform(convert_bam_to_fasta,
           suffix('.fasta'),
           add_inputs([backtranslate_alignment]),
           '.gapped.fasta')
@must_work()
def insert_gaps_task(infiles, outfile):
    infile, (backtranslated,) = infiles
    ref, *seqs = list(SeqIO.parse(infile, 'fasta'))
    ref_gapped = list(r for r in SeqIO.parse(backtranslated, 'fasta')
                      if r.id == ref.id)[0]
    seqs_gapped = (new_record_seq_str(r, insert_gaps(str(ref_gapped.seq),
                                                     str(r.seq),
                                                     '-', '-', skip=1))
                   for r in seqs)
    SeqIO.write(seqs_gapped, outfile, "fasta")


@jobs_limit(n_local_jobs, local_job_limiter)
@active_if(globals_.config.getboolean('Tasks', 'align_full'))
@merge(insert_gaps_task, 'all_timepoints.aligned.fasta')
@must_work(seq_ids=True)
def merge_all_timepoints(infiles, outfile):
    cat_files(infiles, outfile)


def make_alignment_pipeline(name=None):
    if name is None:
        name = "alignment_pipeline"
    pipeline = Pipeline(name)

    n_local_jobs, n_remote_jobs = n_jobs()

    filter_fastq_task = pipeline.transform(filter_fastq,
                                       start_files, suffix(".fastq"), add_inputs([write_config]), '.qfilter.fasta')
    filter_fastq_task.jobs_limit(n_local_jobs, local_job_limiter)

    filter_contaminants_task = pipeline.transform(filter_fastq,
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
                                      output='{path[0]}/{basename[0]}.clusters/cluster_*.raw.fasta')
    cluster_task.jobs_limit(n_remote_jobs, remote_job_limiter)
    cluster_task.mkdir(sort_by_length, suffix('.fasta'), '.clusters')


    select_clusters_task = pipeline.transform(select_clusters,
                                              input=cluster_task,
                                              filter=suffix('.fasta'),
                                              output='.keep.fasta')
    select_clusters_task.jobs_limit(n_local_jobs, local_job_limiter)


    align_clusters_task = pipeline.transform(align_clusters,
                                             input=select_clusters_task,
                                             filter=suffix('.fasta'),
                                             output='.aligned.fasta')
    align_clusters_task.jobs_limit(n_remote_jobs, remote_job_limiter)


    cluster_consensus_task = pipeline.transform(cluster_consensus,
                                                input=align_clusters_task,
                                                filter=suffix('.fasta'),
                                                output='.cons.fasta')
    cluster_consensus_task.jobs_limit(n_local_jobs, local_job_limiter)


    cat_clusters_task = pipeline.collate(cat_clusters,
                                         input=cluster_consensus_task,
                                         filter=formatter(),
                                         output='{subdir[0][0]}.cons.fasta')
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


    add_copynumber_task = pipeline.collate(add_copynumber,
                                           input=[perfect_orfs_task, make_individual_dbs_task, shift_correction_task],
                                           filter=formatter(r'(?P<NAME>.+).qfilter'),
                                           output='{NAME[0]}.copynumber.fasta',
                                           extra='{NAME[0]}')
    add_copynumber_task.jobs_limit(n_remote_jobs, remote_job_limiter)


    cat_all_perfect_task = pipeline.merge(cat_all_perfect,
                                          input=add_copynumber_task,
                                          output="all_perfect_orfs.fasta")
    cat_all_perfect_task.jobs_limit(n_local_jobs, local_job_limiter)


    translate_perfect_task = pipeline.transform(translate_perfect,
                                                input=cat_all_perfect_task,
                                                filter=suffix('.fasta'),
                                                output='.translated.fasta')
    translate_perfect_task.jobs_limit(n_local_jobs, local_job_limiter)


    # output was hyphy_data/input/merged.prot
    # FIXME: copy or move it there in hyphy sub-pipeline
    codon_align_perfect_task = pipeline.transform(codon_align_perfect,
                                                  input=translate_perfect_task,
                                                  filter=formatter(),
                                                  output='merged.prot')
    codon_align_perfect_task.jobs_limit(n_local_jobs, local_job_limiter)
    codon_align_perfect_task.posttask(partial(pause, 'merged.prot'))


    backtranslate_alignment_task = pipeline.merge(backtranslate_alignment,
                                              input=[cat_all_perfect_task,
                                                     codon_align_perfect_task],
                                              output='all_backtranslated.fas')
    backtranslate_alignment_task.jobs_limit(n_local_jobs, local_job_limiter)


    degap_backtranslated_alignment_task = pipeline.transform(degap_backtranslated_alignment,
                                                             input=backtranslate_alignment_task,
                                                             filter=formatter(),
                                                             output='all_perfect_orfs_degapped.fasta')
    degap_backtranslated_alignment_task.jobs_limit(n_local_jobs, local_job_limiter)


    if globals_.config.getboolean('Tasks', 'align_full'):
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
                                                output='{NAME[0]}.alignments/combined.*.unaligned.fasta',
                                                extra='{NAME[0]}')
        combine_pairs_task.jobs_limit(n_local_jobs, local_job_limiter)
        combine_pairs_task.mkdir(full_timestep_pairs, suffix('.pairs.txt'), '.alignments')

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


        insert_gaps_task = pipeline.transform(insert_gaps,
                                              input=convert_bam_to_fasta_task,
                                              filter=suffix('.fasta'),
                                              add_inputs=add_inputs(backtranslate_alignment),
                                              output='.gapped.fasta')
        insert_gaps_task_task.jobs_limit(n_local_jobs, local_job_limiter)


        merge_all_timepoints_task = pipeline.merge(merge_all_timepoints,
                                                   input=insert_gaps_task,
                                                   output='all_timepoints.aligned.fasta')
        merge_all_timepoints_task.jobs_limit(n_local_jobs, local_job_limiter)

    return pipeline
