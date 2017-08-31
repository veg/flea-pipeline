#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
              FLEA (Full-Length Envelope Analyzer) Pipeline
========================================================================================
 FLEA Pipeline. Started March 2016.
 #### Homepage / Documentation
 https://github.com/veg/flea-pipeline
----------------------------------------------------------------------------------------
*/

// TODO: parallelize filter_fastx and other scripts, if it speeds them up
// TODO: parallelize hyphy scripts, if possible
// TODO: fix MDS parallelism
// TODO: set cpus, threads, and time for each task

// TODO: allow starting from aligned inputs
// TODO: nextflow repo format
// TODO: add file extensions
// TODO: Singularity container


params.infile = "$HOME/flea/data/P018/data/metadata"
params.results_dir = "results"


// TODO: how to avoid duplicating?
Channel.fromPath(params.infile)
    .into { metadata_1; metadata_2; metadata_3; metadata_4; metadata_5 }


// read input metadata into tuples
input_files = []
infile = file(params.infile)

infile
    .readLines()
    .each {
        (filename, timepoint_label, date) = it.trim().split()
        tpfile = file(infile.parent / filename)
        mytuple = tuple(tpfile, timepoint_label)
        input_files.add(mytuple)
    }


input_channel = Channel.from(input_files)

/* ************************************************************************** */
/* QUALITY SUB-PIPELINE */

// compute final min/max QCS length from reference length
Channel.fromPath( params.reference_dna )
    .splitFasta( record: [seqString: true ] )
    .map { record -> record.seqString.length() }
    .take(1)
    .into { reflen_1; reflen_2 }

min_qcs_len = reflen_1
    .map { n -> Math.round(n * params.qcs_length_coeff.toBigDecimal()) }

max_qcs_len = reflen_2
    .map { n -> Math.round(n * (2.0 - params.qcs_length_coeff.toBigDecimal())) }


process quality_pipeline {

    cpus params.threads
    time params.slow_time
    publishDir params.results_dir

    input:
    set 'ccs.fastq', label from input_channel
    each minlen from min_qcs_len
    each maxlen from max_qcs_len

    output:
    set '*qcs.fastq', label into qcs_final_1, qcs_final_2

    """
    # filter by quality score
    ${params.usearch} -fastq_filter ccs.fastq \
      -fastq_maxee_rate ${params.max_error_rate} \
      -threads ${params.threads} \
      -fastq_qmax ${params.qmax} -fastq_minlen ${minlen} \
      -fastqout qfiltered \
      -relabel "${label}_ccs_"

    # trim ends
    ${params.python} ${params.script_dir}/trim.py \
      --jobs ${params.threads} --fastq qfiltered trimmed

    # filter runs
    ${params.python} ${params.script_dir}/filter_fastx.py \
      runs fastq fastq ${params.run_length} < trimmed > no_runs

    # filter against contaminants
    ${params.usearch} -usearch_global no_runs \
      -qmask none \
      -db ${params.contaminants_db} -id ${params.contaminant_identity} -strand both \
      -threads ${params.threads} -notmatchedfq uncontam

    # filter against reference
    ${params.usearch} -usearch_global uncontam -db ${params.reference_db} \
      -qmask none \
      -id ${params.reference_identity} \
      -userfields qstrand+tstrand+caln \
      -top_hit_only -strand both \
      -maxaccepts ${params.max_accepts} -maxrejects ${params.max_rejects} \
      -threads ${params.threads} \
      -matchedfq matches_db -userout userout

    # trim terminal gaps
    ${params.python} ${params.script_dir}/trim_terminal_gaps.py matches_db userout filtered

    # length filter
    ${params.python} ${params.script_dir}/filter_fastx.py \
      length fastq fastq \
      ${minlen} ${maxlen} \
      < filtered > ${label}.qcs.fastq
    """
}

/* ************************************************************************** */
/* CONSENSUS SUB-PIPELINE */

process consensus_pipeline {

    cpus params.threads
    time params.crazy_time

    input:
    set 'qcs', label from qcs_final_1

    output:
    file filtered_hqcs into hqcs_files
    file cnfile into cnfiles

    shell:
    '''
    # cluster
    !{params.usearch} -cluster_fast qcs -id !{params.cluster_identity} \
      -qmask none \
      -sort length \
      -maxaccepts !{params.max_accepts} -maxrejects !{params.max_rejects} \
      -top_hit_only -minsl !{params.min_length_ratio} \
      -threads !{params.threads} \
      -clusters raw_cluster_

    # function sample clusters and do mafft consensus
    doconsensus() {
        !{params.python} !{params.script_dir}/filter_fastx.py \
          sample fasta fasta \
          !{params.min_cluster_size} !{params.max_cluster_size} \
          < $1 > ${1}.sampled

        if [ -s ${1}.sampled ]
        then
            !{params.mafft} --ep 0.5 --quiet --preservecase \
            ${1} > ${1}.aligned

             number=$(echo $1 | cut -d '_' -f 3)

             !{params.python} !{params.script_dir}/DNAcons.py -o ${1}.consensus \
                --id !{label}_consensus_${number} ${1}.aligned
        fi

    }
    export -f doconsensus

    # run in parallel
    !{params.parallel} -j !{params.threads} 'doconsensus {}' ::: raw_cluster_*

    cat *.consensus > all_consensus

    ######
    # shift correction

    # compute inframe
    !{params.python} !{params.script_dir}/filter_fastx.py \
      inframe fasta fasta false < all_consensus > inframe

    # unique seqs
    !{params.usearch} -fastx_uniques inframe -fastaout uniques \
      -threads !{params.threads}

    # search
    !{params.usearch} -usearch_global all_consensus -db uniques \
      -qmask none \
      -id !{params.reference_identity} \
      -top_hit_only -strand both \
      -maxaccepts !{params.max_accepts} -maxrejects !{params.max_rejects} \
      -threads !{params.threads} \
      -fastapairs pairfile -userout calnfile -userfields caln

    # shift correction
    !{params.python} !{params.script_dir}/correct_shifts.py \
      --del-strategy=reference \
      --calns=calnfile \
      pairfile shift_corrected

    # compute inframe
    !{params.python} !{params.script_dir}/filter_fastx.py \
      inframe fasta fasta true < shift_corrected > inframe

    # unique seqs
    !{params.usearch} -fastx_uniques inframe -fastaout shifted_uniques \
      -threads !{params.threads}

    ########
    # compute copynumbers

    # usearch for pairs
    !{params.usearch} -usearch_global qcs -db shifted_uniques \
      -qmask none \
      -userout pairfile -userfields query+target \
      -id !{params.copynumber_identity} \
      -top_hit_only -strand both \
      -maxaccepts !{params.max_accepts} -maxrejects !{params.max_rejects} \
      -maxqt !{params.cn_max_length_ratio} \
      -threads !{params.threads}

    # write copynumber file
    !{params.python} !{params.script_dir}/write_copynumbers.py \
      < pairfile > cnfile

    # filter out HQCS with 0 copynumber
    !{params.python} !{params.script_dir}/filter_fastx.py \
      copynumber fasta fasta cnfile \
      < shifted_uniques > filtered_hqcs
    '''
}

process merge_timepoints {

    executor 'local'

    publishDir params.results_dir

    input:
    file 'hqcs*' from hqcs_files.collect()
    file 'cn*' from cnfiles.collect()

    output:
    file 'hqcs.fasta' into merged_hqcs_out

    """
    cat hqcs* > merged_hqcs
    cat cn* > merged_copynumbers

    # add copynumbers to ids, for evo_history
    ${params.python} ${params.script_dir}/add_copynumber_to_id.py \
      merged_hqcs merged_copynumbers hqcs.fasta
    """
}

/* ************************************************************************** */
/* ALIGNMENT SUB-PIPELINE */

process alignment_pipeline {

    time params.slow_time
    cpus params.threads
    time params.slow_time

    publishDir params.results_dir

    input:
    file 'hqcs' from merged_hqcs_out

    output:
    file 'msa.fasta' into msa_out
    file 'msa.aa.fasta' into msa_aa_out

    """
    ${params.python} ${params.script_dir}/translate.py \
      < hqcs > hqcs_protein

    ${params.mafft} --ep 0.5 --quiet --preservecase \
      --thread ${params.threads} \
      hqcs_protein > msa.aa.fasta

    ${params.python} ${params.script_dir}/backtranslate.py \
      msa.aa.fasta hqcs msa.fasta
    """
}

/* ************************************************************************** */
/* ANALYSIS SUB-PIPELINE */

process dates_json_task {

    executor 'local'

    publishDir params.results_dir

    input:
    file 'metadata' from metadata_1

    output:
    file 'dates.json' into dates_json_out

    """
    #!${params.python}
    import json
    from flea.util import get_date_dict

    d = get_date_dict('metadata')
    result = dict((v, k) for k, v in d.items())
    with open('dates.json', 'w') as handle:
        json.dump(result, handle, separators=(",\\n", ":"))
    """
}

process copynumbers_json {

    time params.crazy_time
    publishDir params.results_dir

    input:
    file 'msa' from msa_out

    output:
    file 'copynumbers.json' into copynumbers_json_out

    """
    #!${params.python}

    import json
    from Bio import SeqIO
    from flea.util import id_to_cn

    records = SeqIO.parse('msa', 'fasta')
    d = dict((r.id, id_to_cn(r.id)) for r in records)
    with open('copynumbers.json', 'w') as handle:
        json.dump(d, handle, separators=(",\\n", ":"))
    """
}

process get_oldest_label {

    executor 'local'

    input:
    file 'metadata' from metadata_2

    output:
    stdout oldest_label

    """
    #!${params.python}

    import sys
    from flea.util import get_date_dict

    d = get_date_dict('metadata')
    sys.stdout.write(sorted(d, key=d.get)[0])
    """
}

process mrca {

    time params.slow_time

    input:
    file 'msa' from msa_out
    val oldest_label

    output:
    file mrca into mrca_1, mrca_2, mrca_3
    file mrca_translated into mrca_translated_1, mrca_translated_2

    """
    ${params.python} ${params.script_dir}/filter_fastx.py \
      prefix fasta fasta ${oldest_label} \
      < msa > oldest_seqs

    ${params.python} ${params.script_dir}/DNAcons.py \
      --keep-gaps --codon --id MRCA \
      -o mrca \
      --copynumbers \
      oldest_seqs

    ${params.python} ${params.script_dir}/translate.py --gapped \
      < mrca > mrca_translated
    """
}

// TODO: why do we have to duplicate outputs here?
process add_mrca {

    executor 'local'

    input:
    file 'mrca' from mrca_1
    file 'msa' from msa_out

    output:
    file 'msa_with_mrca' into msa_with_mrca_1, msa_with_mrca_2

    "cat mrca msa > msa_with_mrca"
}

process fasttree {

    time params.slow_time

    cpus params.threads
    time params.slow_time

    input:
    file 'msa' from msa_with_mrca_1

    output:
    file tree

    "OMP_NUM_THREADS=${params.threads} ${params.fasttree} -gtr -nt msa > tree"
}

process reroot {

    time params.fast_time

    input:
    file tree

    output:
    file rooted_tree into rooted_tree_1, rooted_tree_2

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

process tree_json {

    executor 'local'

    publishDir params.results_dir

    input:
    file 'tree' from rooted_tree_1

    output:
    file 'trees.json' into tree_json_out

    """
    #!${params.python}
    import json

    with open('tree') as handle:
        newick_string = handle.read()
    result = {'tree': newick_string}
    with open('trees.json', 'w') as handle:
        json.dump(result, handle, separators=(",\\n", ":"))
    """
}

process js_divergence {

    time params.slow_time

    publishDir params.results_dir

    input:
    file 'msa' from msa_aa_out
    file 'metadata' from metadata_3

    output:
    file 'js_divergence.json' into js_divergence_json

    """
    ${params.python} ${params.script_dir}/js_divergence.py \
      msa metadata js_divergence.json
    """
}

process manifold_embedding {

    time params.crazy_time

    publishDir params.results_dir

    input:
    file 'msa' from msa_out

    output:
    file 'manifold.json' into manifold_json_out

    """
    ${params.usearch} -calc_distmx msa -distmxout dmatrix \
      -format tabbed_pairs

    ${params.python} ${params.script_dir}/manifold_embed.py \
      --n-jobs 1 dmatrix manifold.json
    """
}

// TODO: avoid full paths
// TODO: why do command line arguments not work here?
process reconstruct_ancestors {

    time params.slow_time

    input:
    file 'msa' from msa_with_mrca_2
    file 'msa_aa' from msa_aa_out
    file 'tree' from rooted_tree_2

    output:
    file msa_aa_ancestors into msa_aa_ancestors_out

    shell:
    '''
    echo $(pwd)/msa >> stdin
    echo $(pwd)/tree >> stdin
    echo $(pwd)/ancestors >> stdin
    echo HKY85 >> stdin
    echo 2 >> stdin
    !{params.hyphy} !{params.hyphy_dir}/reconstructAncestors.bf < stdin

    !{params.python} !{params.script_dir}/translate.py --gapped \
      < ancestors > translated

    cat msa_aa translated > msa_aa_ancestors
    '''
}

process coordinates_json {

    time params.fast_time

    cpus params.threads
    time params.slow_time

    publishDir params.results_dir

    input:
    file 'mrca' from mrca_translated_1

    output:
    file 'coordinates.json' into coordinates_json_out

    """
    cat mrca ${params.reference_protein} > pair.fasta

    ${params.mafft} --ep 0.5 --quiet --preservecase \
      --thread ${params.threads} \
      pair.fasta > aligned.fasta

    ${params.python} ${params.script_dir}/coordinates_json.py \
      mrca aligned.fasta ${params.reference_coordinates} coordinates.json
    """
}

process sequences_json {

    time params.fast_time

    publishDir params.results_dir

    input:
    file 'msa' from msa_aa_ancestors_out
    file 'mrca' from mrca_translated_2
    file 'coordinates.json' from coordinates_json_out
    file 'metadata' from metadata_4

    output:
    file 'sequences.json' into sequences_json_out

    """
    ${params.python} ${params.script_dir}/sequences_json.py \
      msa mrca coordinates.json metadata \
      ${params.reference_protein} ${params.reference_coordinates} \
      sequences.json
    """
}

process replace_stop_codons {

    time params.fast_time

    input:
    file 'msa' from msa_out

    output:
    file no_stops into msa_no_stops

    """
    ${params.python} ${params.script_dir}/filter_fastx.py \
      stop_codons fasta fasta < msa > no_stops
    """
}

// TODO: why do we need to split output here, but not elsewhere?
process seq_dates {

    time params.fast_time

    input:
    file 'msa' from msa_out
    file 'metadata' from metadata_5

    output:
    file dates into seq_dates_1, seq_dates_2

    """
    #!${params.python}
    import json
    from Bio import SeqIO
    from flea.util import get_date_dict

    date_dict = get_date_dict('metadata')
    records = list(SeqIO.parse('msa', "fasta"))
    with open('dates', "w") as handle:
        outdict = {r.id: date_dict[r.id.split('_')[0]] for r in records}
        json.dump(outdict, handle, separators=(",\\n", ":"))
    """
}

process region_coords {

    time params.fast_time

    publishDir params.results_dir

    input:
    file 'mrca' from mrca_2

    output:
    file 'region_coords.json' into region_coords_json

    shell:
    '''
    !{params.hyphy} !{params.hyphy_dir}/HXB2partsSplitter.bf \
      $(pwd)/mrca $(pwd)/region_coords.json
    '''
}


process evo_history {

    time params.crazy_time

    input:
    file 'no_stops' from msa_no_stops
    file 'dates' from seq_dates_1
    file 'region_coords' from region_coords_json

    output:
    file 'rates_pheno.tsv' into rates_pheno

    shell:
    '''
    !{params.hyphy} !{params.hyphy_dir}/obtainEvolutionaryHistory.bf \
      $(pwd)/no_stops $(pwd)/dates $(pwd)/region_coords $(pwd)/rates_pheno.tsv
    '''
}

process rates_pheno_json {

    time params.fast_time
    publishDir params.results_dir

    input:
    file rates_pheno

    output:
    file 'rates_pheno.json' into rates_pheno_json_out

    """
    #!${params.python}

    import csv
    import json

    with open('rates_pheno.tsv') as h:
        lines = list(csv.reader(h, delimiter='\t'))
    result = {}
    keys, rest = lines[0], lines[1:]
    result = list(dict((key, value) for key, value in zip(keys, line))
                  for line in rest)
    with open('rates_pheno.json', 'w') as h:
        json.dump(result, h, separators=(",\\n", ":"))
    """
}

process fubar {

    time params.crazy_time
    publishDir params.results_dir

    input:
    file 'no_stops' from msa_no_stops
    file 'dates' from seq_dates_2
    file 'mrca' from mrca_3

    output:
    file 'rates.json' into rates_json

    shell:
    '''
    !{params.hyphy} !{params.hyphy_dir}/runFUBAR.bf \
      $(pwd)/no_stops $(pwd)/dates $(pwd)/mrca $(pwd) $(pwd)/rates.json
    '''
}

/* ************************************************************************** */
/* DIAGNOSIS SUB-PIPELINE */

process diagnose {
    when:
    params.do_diagnosis

    cpus params.threads
    time params.crazy_time

    publishDir params.results_dir

    input:
    file 'qcs*' from qcs_final_2.map{ it[0] }.collect()
    file 'hqcs' from merged_hqcs_out
    file 'hqcs_msa' from msa_out
    file 'hqcs_aa_msa' from msa_aa_out

    output:
    file 'diagnosis_results/*' into diagnosis_results

    shell:
    '''
    # cat qcs
    cat qcs* > all_qcs

    # convert to fasta
    !{params.python} !{params.script_dir}/filter_fastx.py \
      convert fastq fasta < all_qcs > qcs.fasta

    # usearch for pairs
    !{params.usearch} -usearch_global qcs.fasta -db hqcs \
      -qmask none \
      -userout pairfile -userfields query+target \
      -id !{params.qcs_to_hqcs_identity} \
      -top_hit_only -strand both \
      -maxaccepts !{params.max_accepts} -maxrejects !{params.max_rejects} \
      -maxqt !{params.max_qcs_length_ratio} \
      -threads !{params.threads}

    # combine pairs
    mkdir alignments
    !{params.python} !{params.script_dir}/pairs_to_fasta.py \
      qcs.fasta hqcs pairfile alignments

    # align and insert gaps
    !{params.parallel} -j !{params.threads} \
      '!{params.bealign} -a codon {} {}.bam' ::: alignments/*unaligned.fasta

    !{params.parallel} -j !{params.threads} \
      '!{params.bam2msa} {} {}.fasta' ::: alignments/*.bam

    !{params.parallel} -j !{params.threads} \
      '!{params.python} !{params.script_dir}/insert_gaps.py {} hqcs_msa {}.gapped' ::: alignments/*.bam.fasta

    # merge all
    cat alignments/*gapped > qcs_msa

    # replace gapped codons
    !{params.python} !{params.script_dir}/filter_fastx.py \
      gap_codons fasta fasta < qcs_msa > qcs_msa_nogaps

    !{params.python} !{params.script_dir}/translate.py --gapped \
      < qcs_msa_nogaps > qcs_msa_aa

    # run diagnosis
    mkdir diagnosis_results
    !{params.python} !{params.script_dir}/diagnose.py \
      hqcs_aa_msa qcs_msa_aa diagnosis_results
    '''
}