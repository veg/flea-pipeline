#!/usr/bin/env nextflow

/*
 * The FLEA pipeline
 *
 */

// TODO: run on TORQUE
// TODO: write consensus, alignment, analysis, diagnosis pipelines
// TODO: allow starting from aligned inputs
// TODO: do not hardcode min/max lengths
// TODO: sub-pipelines
// TODO: nextflow repo format

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

    """${params.python} ${params.script_dir}/trim.py \
      --jobs ${params.threads} --fastq qfiltered trimmed
    """
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

filter_length_outs.into { qcs_outs_for_cluster; qcs_outs_for_copynumber }

process cluster {
    input:
    set 'lenfiltered', label from qcs_outs_for_cluster

    output:
    set 'raw_cluster_*', label into cluster_outs

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

/*
 * Shift-corrected, inframe, unique HQCS sequences.
 */
process shift_correction {
    input:
    file all_consensus

    output:
    file shifted_uniques

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
      -maxaccepts ${params.max_accepts} -maxrejects ${params.max_rejects} \
      -threads ${params.threads} \
      -fastapairs pairfile -userout calnfile -userfields caln

    # shift correction
    ${params.python} ${params.script_dir}/correct_shifts.py \
      --del-strategy=reference \
      --calns=calnfile \
      pairfile shift_corrected

    # compute inframe
    ${params.python} ${params.script_dir}/filter_fastx.py \
      inframe fasta fasta true < shift_corrected > inframe

    # unique seqs
    ${params.usearch} -fastx_uniques inframe -fastaout shifted_uniques \
      -threads ${params.threads}
      
    """
}

// get hqcs and qcs sequences for this time point
shifted_uniques
    .merge( qcs_outs_for_copynumber.map { it[0] } ) { hqcs, ccs -> [hqcs, ccs] }
    .set { hqcs_qcs_pairs }

/*
 * Compute pairs for copynumbers.
 */
process copynumbers {
    input:
    set 'hqcs', 'qcs' from hqcs_qcs_pairs

    output:
    set 'final_hqcs', 'cnfile' into hqcs_copynumber_pairs

    """
    # usearch for pairs
    ${params.usearch} -usearch_global qcs -db hqcs \
      -userout pairfile -userfields query+target \
      -id ${params.copynumber_identity} \
      -top_hit_only -strand both \
      -maxaccepts ${params.max_accepts} -maxrejects ${params.max_rejects} \
      -maxqt ${params.cn_max_length_ratio} \
      -threads ${params.threads}

    # write copynumber file
    ${params.python} ${params.script_dir}/write_copynumbers.py \
      < pairfile > cnfile

    # filter out HQCS with 0 copynumber
    ${params.python} ${params.script_dir}/filter_fastx.py \
      copynumber fasta fasta cnfile \
      < hqcs > final_hqcs
    """
}

hqcs_files = Channel.create()
cn_files = Channel.create()

hqcs_copynumber_pairs
    .separate( hqcs_files, cn_files ) { x -> [x[0], x[1]] }

process merge_hqcs {
    input:
    file 'hqcs*' from hqcs_files.collect()

    output:
    file all_hqcs

    """
    cat hqcs* > all_hqcs
    """
}

process merge_copynumbers {
    input:
    file 'cn*' from cn_files.collect()

    output:
    file all_cns

    """
    cat cn* > all_cns
    """
}

/* ************************************************************************** */
/* ALIGNMENT SUB-PIPELINE */

all_hqcs.into { hqcs_for_translate; hqcs_for_backtranslate }

process translate_hqcs {
    input:
    file hqcs_for_translate

    output:
    file hqcs_translated

    """
    ${params.python} ${params.script_dir}/translate.py \
      < all_hqcs > hqcs_translated
    """
}

process align_hqcs {
    input:
    file hqcs_translated

    output:
    file hqcs_aligned

    """
    ${params.mafft} --ep 0.5 --quiet --preservecase \
      hqcs_translated > hqcs_aligned
    """
}

hqcs_for_backtranslate
    .merge ( hqcs_aligned ) { dna, aa -> [dna, aa] }
    .set { dna_aa_pairs }

// TODO: threaded backtranslate
process backtranslate_hqcs {
    input:
    set 'dna', 'aa' from dna_aa_pairs

    output:
    file backtranslated

    """
    ${params.python} ${params.script_dir}/backtranslate.py \
      aa dna backtranslated
    """
}

// FIXME: why does the pipeline hang here???

/* ************************************************************************** */
/* ANALYSIS SUB-PIPELINE */

// TODO: combine into fewest possible tasks
// TODO: thread as many as possible
// TODO: all .json files
// TODO: .zip file

// duplicate all the inputs needed by this part of the pipeline

/*
 * evo history needs copynumbers in the names
 */
process add_cn_to_ids {
    input:
    set 'msa', 'cns' from XXX

    output:
    msa_with_cn

    """
    #!${params.python}
    from Bio import SeqIO
    from flea.util import parse_copynumbers, replace_id
    
    id_to_cn = parse_copynumbers('cns')
    records = SeqIO.parse('msa', 'fasta')
    result = (replace_id(r, id_with_cn(r.id, id_to_cn[r.id]))
              for r in records)
    SeqIO.write(result, 'msa_with_cn', "fasta")
    """
}

process mrca {
    input:

    output:

    """
    #!${params.python}

    alignment_file, copynumber_file = infiles
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(globals_.timepoints, key=strptime)
    oldest_records_filename = os.path.join(pipeline_dir, 'oldest_sequences.fasta')

    records = SeqIO.parse(alignment_file, "fasta")
    oldest_records = (r for r in records if r.id.startswith(oldest_timepoint.label))
    SeqIO.write(oldest_records, oldest_records_filename, "fasta")

    return mrca(oldest_records_filename, copynumber_file, outfile)
    """
}

process add_mrca {

}

process fasttree {

}

process reroot {

}

process translate_mrca {

}

process translate_msa {

}

process js_divergence {

}

process distance_matrix {

}

process manifold_embedding {

}

process reconstruct_ancestors {

}

process translate_ancestors {

}

process replace_stop_codons {

}

process dates_for_evo_history {

}

process evo_history {

}

process fubar {

}


/* ************************************************************************** */
/* DIAGNOSIS SUB-PIPELINE */

// TODO: make this optional
