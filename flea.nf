#!/usr/bin/env nextflow

/*
 * The FLEA pipeline
 *
 */

// TODO: get it to run on multiple inputs in parallel
// TODO: run on TORQUE
// TODO: write consensus, alignment, analysis, diagnosis pipelines
// TODO: allow starting from aligned inputs
// TODO: do not hardcode min/max lengths

params.infile = "$HOME/flea/data/P018/data/metadata"

// read input metadata into tuples
input_files = []
tp_to_date = {}
infile = file(params.infile)

infile
    .readLines()
    .each {
        (filename, timepoint_label, date) = it.trim().split()
        tpfile = file(infile.parent / filename)
        mytuple = tuple(tpfile, timepoint_label)
        input_files.add(mytuple)
        tp_to_date[timepoint_label] = date
    }

input_channel = Channel.from(input_files)

/* ************************************************************************** */
/* QUALITY SUB-PIPELINE */

process quality_filter {
    input:
    set 'ccs.fastq', label from input_channel

    output:
    set 'qfiltered', label into quality_outs

    """
    ${params.usearch} -fastq_filter ccs.fastq \
      -fastq_maxee_rate ${params.max_error_rate} \
      -threads ${params.threads} \
      -fastq_qmax ${params.qmax} -fastq_minlen ${params.min_length} \
      -fastqout qfiltered \
      -relabel "${label}_ccs_"
    """
}

process trim_ends {
    input:
    set 'qfiltered', label from quality_outs

    output:
    set 'trimmed', label into trim_ends_outs

    "${params.python} ${params.script_dir}/trim.py --jobs ${params.threads} --fastq qfiltered trimmed"
}

process filter_runs {
    input:
    set 'trimmed', label from trim_ends_outs

    output:
    set 'no_runs', label into filter_runs_outs

    """
    ${params.python} ${params.script_dir}/filter_fastx.py \
      runs fastq fastq ${params.run_length} < trimmed > no_runs
    """
}

process filter_db {
    input:
    set 'no_runs', label from filter_runs_outs

    output:
    set 'filtered', label into filter_db_outs

    """
    ${params.usearch} -usearch_global no_runs \
      -db ${params.contaminants_db} -id ${params.contaminant_identity} -strand both \
      -threads ${params.threads} -notmatchedfq uncontam

    ${params.usearch} -usearch_global uncontam -db ${params.reference_db} \
      -id ${params.reference_identity} \
      -userfields qstrand+tstrand+caln \
      -top_hit_only -strand both \
      -maxaccepts ${params.max_accepts} -maxrejects ${params.max_rejects} \
      -threads ${params.threads} \
      -matchedfq matches_db -userout userout

     ${params.python} ${params.script_dir}/trim_terminal_gaps.py matches_db userout filtered
     """
}

// FIXME: do not hardcode min/max length
process filter_length {
    input:
    set 'filtered', label from filter_db_outs

    output:
    set 'lenfiltered', label into filter_length_outs

    """
    ${params.python} ${params.script_dir}/filter_fastx.py \
      length fastq fastq \
      ${params.min_length} ${params.max_length} \
      < filtered > lenfiltered
    """
}

/* ************************************************************************** */
/* CONSENSUS SUB-PIPELINE */

process cluster {
    input:
    set 'lenfiltered', label from filter_length_outs

    output:
    set "raw_cluster_*", label into cluster_outs

    """
    ${params.usearch} -cluster_fast lenfiltered -id ${params.cluster_identity} \
      -sort length \
      -maxaccepts ${params.max_accepts} -maxrejects ${params.max_rejects} \
      -top_hit_only -minsl ${params.min_length_ratio} \
      -threads ${params.threads} \
      -clusters raw_cluster_
    """
}

process consensus {
    input:
    set 'raw*', label from cluster_outs

    output:
    file all_consensus

    shell:
    '''
    for f in raw*; do
        !{params.python} !{params.script_dir}/filter_fastx.py \
          sample fasta fasta \
          !{params.min_cluster_size} !{params.max_cluster_size} \
          < $f > ${f}.sampled

        if [ -s ${f}.sampled ]
	then
            !{params.mafft} --ep 0.5 --quiet --preservecase \
            ${f} > ${f}.aligned

             number=$(echo $f | cut -d '_' -f 2)

             !{params.python} !{params.script_dir}/DNAcons.py -o ${f}.consensus \
                --id !{label}_consensus_${number} ${f}.aligned
        fi
    done
    cat *.consensus > all_consensus
    '''
}

process shift_correction {
    input:
    file all_consensus

    output:
    file shift_corrected

    """
    # compute inframe
    ${params.python} ${params.script_dir}/filter_fastx.py \
      inframe fasta fasta false < all_consensus > inframe

    # unique seqs
    ${params.usearch} -fastx_uniques inframe -fastaout uniques \
      -threads ${params.threads}

    # search
    ${params.usearch} -usearch_global all_consensus -db uniques \
      -id ${params.reference_identity} \
      -top_hit_only -strand both \
      -maxaccepts {params.max_accepts} -maxrejects {params.max_rejects} \
      -threads ${params.threads} \
      -fastapairs pairfile -userout calnfile -userfields caln

    # shift correction
    ${params.python} ${params.script_dir}/correct_shifts.py \
      --del-strategy=reference \
      --calns=calnfile \
      pairfile shift_corrected
    """
}

// TODO: duplicate qcs file after quality sub-pipeline, so it can be used here.
/*
 * Get inframe, unique HQCS sequences.
 * Compute pairs for copynumbers.
 */
process finalize_hqcs {
    input:
    file shift_corrected

    output:
    file final_hqcs

    """
    # compute inframe
    ${params.python} ${params.script_dir}/filter_fastx.py \
      inframe fasta fasta true < shift_corrected > inframe

    # unique seqs
    ${params.usearch} -fastx_uniques inframe -fastaout uniques \
      -threads ${params.threads}

    # usearch for pairs
    ${params.usearch} -usearch_global qcsfile -db uniques \
      -userout pairfile -userfields query+target \
      -id ${params.copynumber_identity} \
      -top_hit_only -strand both \
      -maxaccepts {params.max_accepts} -maxrejects {params.max_rejects}
      -maxqt ${params.cn_max_length_ratio} \
      -threads {ppn}

    """
}

// TODO: combine pairs to make copynumbers
// TODO: only keep HQCSs with nonzero copynumber.
// TODO: merge all files
