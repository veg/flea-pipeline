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
// TODO: infer singleton proc

params.infile = "$HOME/flea/data/P018/data/metadata"

Channel.fromPath(params.infile)
    .into { metadata_1; metadata_2 }


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

// FIXME: do this programmatically
oldest_label = 'V06'

input_channel = Channel.from(input_files)

/* ************************************************************************** */
/* QUALITY SUB-PIPELINE */

process quality_filter {
    input:
    set 'ccs.fastq', label from input_channel

    output:
    set 'qfiltered', label into quality_out

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
    set 'qfiltered', label from quality_out

    output:
    set 'trimmed', label into trim_ends_out

    """${params.python} ${params.script_dir}/trim.py \
      --jobs ${params.threads} --fastq qfiltered trimmed
    """
}

process filter_runs {
    input:
    set 'trimmed', label from trim_ends_out

    output:
    set 'no_runs', label into filter_runs_out

    """
    ${params.python} ${params.script_dir}/filter_fastx.py \
      runs fastq fastq ${params.run_length} < trimmed > no_runs
    """
}

process filter_db {
    input:
    set 'no_runs', label from filter_runs_out

    output:
    set 'filtered', label into filter_db_out

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
    set 'filtered', label from filter_db_out

    output:
    set 'lenfiltered', label into qcs_final_1, qcs_final_2

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
    set 'lenfiltered', label from qcs_final_1

    output:
    set 'raw_cluster_*', label into cluster_out

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
    set 'raw*', label from cluster_out

    output:
    set 'all_consensus', label into consensus_out

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
    set 'all_consensus', label from consensus_out

    output:
    set 'shifted_uniques', label into shift_correction_out

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

/*
 * Compute pairs for copynumbers.
 *
 * We need to be careful to group files by time point.
 *
 */
process copynumbers {
    input:
    set label, 'hqcs', 'qcs' from shift_correction_out.combine(qcs_final_2, by: 1)

    output:
    file filtered_hqcs into hqcs_files
    file cnfile into cnfiles

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
      < hqcs > filtered_hqcs
    """
}


process merge_hqcs {
    input:
    file 'hqcs' from hqcs_files.collect()

    output:
    file merged_hqcs into merged_hqcs_1, merged_hqcs_2

    """
    cat hqcs* > merged_hqcs
    """
}

process merge_copynumbers {
    input:
    file 'cn*' from cnfiles.collect()

    output:
    file all_cns into all_cns_1, all_cns_2

    """
    cat cn* > all_cns
    """
}

/* ************************************************************************** */
/* ALIGNMENT SUB-PIPELINE */

process translate_hqcs {
    input:
    file 'hqcs' from merged_hqcs_1

    output:
    file hqcs_translated

    """
    ${params.python} ${params.script_dir}/translate.py \
      < hqcs > hqcs_translated
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

// TODO: threaded backtranslate
process backtranslate_hqcs {
    input:
    file 'dna' from merged_hqcs_2
    file 'aa' from hqcs_aligned

    output:
    file msa into msa_1, msa_2

    """
    ${params.python} ${params.script_dir}/backtranslate.py \
      aa dna msa
    """
}

/* ************************************************************************** */
/* ANALYSIS SUB-PIPELINE */

// TODO: combine into fewest possible tasks
// TODO: thread as many as possible
// TODO: all .json files
// TODO: .zip file

// TODO: rewrite as filter_fastx with id prefix
process oldest_seqs {
    input:
    file 'msa' from msa_1

    output:
    file oldest_seqs

    """
    #!${params.python}
    import datetime
    from Bio import SeqIO

    records = SeqIO.parse('msa', "fasta")
    oldest_records = (r for r in records if r.id.startswith("${oldest_label}"))
    SeqIO.write(oldest_records, 'oldest_seqs', "fasta")
    """
}

process mrca {
    input:
    file oldest_seqs
    file 'cns' from all_cns_1

    output:
    file 'mrca' into mrca_1, mrca_2

    """
    ${params.python} ${params.script_dir}/DNAcons.py \
      --keep-gaps --codon --id MRCA \
      -o mrca \
      --copynumbers cns \
      oldest_seqs
    """
}

/*
 * evo history needs copynumbers in the names
 */
process add_cn_to_ids {
    input:
    file 'msa' from msa_2
    file 'cns' from all_cns_2

    output:
    file msa_with_cn into msa_id_1, msa_id_2, msa_id_3, msa_id_4

    """
    #!${params.python}
    from Bio import SeqIO
    from flea.util import parse_copynumbers, replace_id

    def id_with_cn(id_, cn):
        return "{}_cn_{}".format(id_, cn)

    id_to_cn = parse_copynumbers('cns')
    records = SeqIO.parse('msa', 'fasta')
    result = (replace_id(r, id_with_cn(r.id, id_to_cn[r.id]))
              for r in records)
    SeqIO.write(result, 'msa_with_cn', "fasta")
    """
}

process add_mrca {
    input:
    file 'mrca' from mrca_1
    file 'msa' from msa_id_1

    output:
    file msa_with_mrca

    "cat mrca msa > msa_with_mrca"
}

process fasttree {
    input:
    file msa_with_mrca

    output:
    file tree

    "${params.fasttree} -gtr -nt msa_with_mrca > tree"
}

process reroot {
    input:
    file tree

    output:
    file rooted_tree

    """
    #!${params.python}
    from Bio import Phylo

    tree = next(Phylo.parse('tree', 'newick'))
    clade = next(tree.find_clades('MRCA'))
    tree.root_with_outgroup(clade)

    # also rename for HyPhy
    for i, node in enumerate(tree.get_nonterminals()):
        node.confidence = None
        if node.name != 'MRCA':
            node.name = "ancestor_{}".format(i)
    Phylo.write([tree], 'rooted_tree', 'newick')
    """
}

process translate_msa_with_cn {
    input:
    file 'msa' from msa_id_2

    output:
    file msa_translated

    """
    ${params.python} ${params.script_dir}/translate.py --gapped \
      < msa > msa_translated
    """
}

process translate_mrca {
    input:
    file 'mrca' from mrca_2

    output:
    file mrca_translated

    """
    ${params.python} ${params.script_dir}/translate.py --gapped \
      < mrca > mrca_translated
    """
}

process js_divergence {
    input:
    file 'msa' from msa_id_2
    file 'metadata' from metadata_1

    output:
    file 'js_divergence.json' into js_divergence_json

    """
    ${params.python} ${params.script_dir}/js_divergence.py \
      msa metadata js_divergence.json
    """
}

// process distance_matrix {
//     input:
//     file 'msa' from msa_cn_4

//     output:
//     file dmatrix

//     """
//     ${params.usearch} -calc_distmx msa -distmxout dmatrix \
//       -format tabbed_pairs
//     """
// }

// process manifold_embedding {
//     input:
//     file dmatrix

//     output:
//     file 'manifold.json' into manifold_json_out

//     """
//     ${params.python} ${params.script_dir}/manifold_embed.py \
//       --n-jobs ${params.threads} --flip \
//       dmatrix manifold.json
//     """
// }

// // process reconstruct_ancestors {

// // }

// // process translate_ancestors {

// // }

// // process replace_stop_codons {

// // }

// // process dates_for_evo_history {

// // }

// // process evo_history {

// // }

// // process fubar {

// // }

// process dates_json_task {
//     input:
//     file 'metadata' from metadata_2

//     output:
//     file 'dates.json' into dates_json_out

//     """
//     #!${params.python}
//     import json
//     from flea.util import get_date_dict

//     result = get_date_dict('metadata')
//     with open('dates.json', 'w') as handle:
//         json.dump(result, handle, separators=(",\\n", ":"))
//     """
// }

// // /* ************************************************************************** */
// // /* DIAGNOSIS SUB-PIPELINE */

// // // TODO: make this optional
