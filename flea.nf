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


// TODO: update tests
// TODO: tune maxaccepts and maxrejects
// TODO: combine all time points for inframe db for shift correction


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

    tag { label }

    publishDir params.results_dir

    cpus params.cpus
    time params.slow_time

    input:
    set 'ccs.fastq', label from input_channel
    each minlen from min_qcs_len
    each maxlen from max_qcs_len

    output:
    set '*qcs.fastq.gz', label into qcs_final_1, qcs_final_2, qcs_final_3, qcs_final_4

    shell:
    '''
    # filter by quality score
    !{params.usearch} --fastq_filter ccs.fastq \
      --fastqout qfiltered.fastq \
      --fastq_maxee_rate !{params.max_error_rate} \
      --fastq_qmax !{params.qmax} \
      --fastq_minlen !{minlen} \
      --relabel "!{label}_ccs_" \
      --threads !{params.cpus}

    # trim ends
    !{params.python} !{params.script_dir}/trim_tails.py \
      --jobs !{params.cpus} --fastq qfiltered.fastq trimmed.fastq

    # filter runs
    !{params.python} !{params.script_dir}/filter_fastx.py \
      runs fastq fasta !{params.run_length} < trimmed.fastq > no_runs.fasta

    # filter against contaminants
    !{params.usearch} --usearch_global no_runs.fasta \
      --db !{params.contaminants_db} \
      --notmatched uncontam.fasta \
      --id !{params.contaminant_identity} \
      --strand both \
      --qmask none \
      --top_hit_only \
      --maxrejects !{params.max_rejects} \
      --threads !{params.cpus}

    # filter against reference
    !{params.usearch} --usearch_global uncontam.fasta \
      --db !{params.reference_db} \
      --userout userout.txt \
      --userfields query+qstrand+tstrand+caln \
      --id !{params.reference_identity} \
      --qmask none \
      --strand both \
      --top_hit_only \
      --maxaccepts !{params.max_accepts} \
      --maxrejects !{params.max_rejects} \
      --threads !{params.cpus}

    # propagate db search to fastq file. trim terminal gaps.
    # use trimmed.fastq since that was the last fastq file in the pipeline.
    # any sequences filtered out of `no_runs` won't make it through the database
    # searches, so it will get filtered out here again.
    !{params.python} !{params.script_dir}/filter_fastx.py \
      userout fastq fastq userout.txt \
      < trimmed.fastq > filtered.fastq

    # length filter
    !{params.python} !{params.script_dir}/filter_fastx.py \
      length fastq fastq \
      !{minlen} !{maxlen} \
      < filtered.fastq > !{label}.qcs.fastq

    # compress all files
    for i in `find . 1 ! -type l | grep -E "\\.fasta$|\\.fastq$|\\.txt$"`; do gzip "$i" ; done
    '''
}

/* ************************************************************************** */
/* CONSENSUS SUB-PIPELINE */


process cluster {

    tag { label }

    publishDir params.results_dir

    cpus params.cpus
    time params.slow_time

    input:
    set 'qcs.fastq.gz', label from qcs_final_1

    output:
    set '*.clusters.uc', label into cluster_out

    shell:
    '''
    zcat qcs.fastq.gz > qcs.fastq

    # cluster
    !{params.usearch} --cluster_fast qcs.fastq \
      -uc !{label}.clusters.uc \
      --id !{params.cluster_identity} \
      --minsl !{params.min_length_ratio} \
      --sort length \
      --top_hit_only \
      --maxaccepts !{params.max_accepts} \
      --maxrejects !{params.max_rejects} \
      --threads !{params.cpus}

    rm -f qcs.fastq
    '''
}


cluster_out
  .phase (qcs_final_2) { it[1] }
  .map { [ it[0][0], it[1][0], it[0][1] ] }
  .set { consensus_input }


process consensus {

    tag { label }

    publishDir params.results_dir

    cpus params.cpus
    time params.slow_time

    input:
    set 'clusters.uc', 'qcs.fastq.gz', label from consensus_input

    output:
    set '*.all_consensus.fasta.gz', label into consensus_out

    shell:
    '''
    zcat qcs.fastq.gz | \
    !{params.python} !{params.script_dir}/cluster_fastq.py \
      --minsize !{params.min_cluster_size} --fasta \
      clusters.uc .

    # function sample clusters and do mafft consensus
    doconsensus() {
        !{params.python} !{params.script_dir}/filter_fastx.py \
          sample fasta fasta \
          !{params.min_cluster_size} !{params.max_cluster_size} \
          < ${1} > ${1}.sampled.fasta

        if [ -s ${1}.sampled.fasta ]
        then
            !{params.mafft} --ep 0.5 --quiet --preservecase \
              ${1}.sampled.fasta > ${1}.sampled.aligned.fasta

	    # get cluster number so we can put it in the record id
            number=$(echo ${1} | cut -d '_' -f 2)

            !{params.python} !{params.script_dir}/DNAcons.py \
	      -o ${1}.consensus.fasta \
              --id !{label}_consensus_${number} ${1}.sampled.aligned.fasta
        fi

    }
    export -f doconsensus

    # run in parallel
    !{params.parallel} -j !{params.cpus} 'doconsensus {}' ::: *_raw.fasta

    cat *.consensus.fasta > !{label}.all_consensus.fasta

    # gzip everything
    rm -f qcs.fastq
    for i in `find . ! -type l | grep -E "\\.fasta$|\\.fastq$|\\.txt$"`; do gzip "$i" ; done
    '''
}

process shift_correction {

    tag { label }

    publishDir params.results_dir

    cpus params.cpus
    time params.slow_time

    input:
    set 'consensus.fasta.gz', label from consensus_out

    output:
    set '*.corrected.inframe.unique.fasta.gz', label into shift_correction_out

    shell:
    '''
    zcat consensus.fasta.gz > consensus.fasta

    # inframe
    !{params.python} !{params.script_dir}/filter_fastx.py \
      inframe fasta fasta false \
      < consensus.fasta > consensus.inframe.fasta

    # unique
    !{params.usearch} --fastx_uniques consensus.inframe.fasta \
      --fastaout consensus.inframe.unique.fasta \
      --threads !{params.cpus}

    # search
    !{params.usearch} --usearch_global consensus.fasta \
      --db consensus.inframe.unique.fasta \
      --fastapairs pairfile.fasta \
      --userout calnfile.txt \
      --userfields caln \
      --top_hit_only \
      --id !{params.reference_identity} \
      --qmask none \
      --strand plus \
      --maxaccepts !{params.max_accepts} \
      --maxrejects !{params.max_rejects} \
      --threads !{params.cpus}

    # shift correction
    !{params.python} !{params.script_dir}/correct_shifts.py \
      --del-strategy=reference \
      --calns=calnfile.txt \
      pairfile.fasta corrected.fasta

    # compute inframe
    !{params.python} !{params.script_dir}/filter_fastx.py \
      inframe fasta fasta true < corrected.fasta > corrected.inframe.fasta

    # unique seqs
    !{params.usearch} --fastx_uniques corrected.inframe.fasta \
      --fastaout !{label}.corrected.inframe.unique.fasta \
      --threads !{params.cpus}

    # gzip everything
    rm -f consensus.fasta
    for i in `find . ! -type l | grep -E "\\.fasta$|\\.fastq$|\\.txt$"`; do gzip -f "$i" ; done
    '''
}

shift_correction_out
  .phase (qcs_final_3) { it[1] }
  .map { [ it[0][0], it[1][0], it[0][1] ] }
  .set { compute_abundances_input }


process compute_abundances {

    tag { label }

    cpus params.cpus
    time params.slow_time

    input:
    set 'hqcs.fasta.gz', 'qcs.fastq.gz', label from compute_abundances_input

    output:
    file 'hqcs.filtered.fasta.gz' into hqcs_files
    file 'abundance_file.txt.gz' into abundance_files

    shell:
    '''
    zcat hqcs.fasta.gz > hqcs.fasta

    # convert to fasta for usearch
    zcat qcs.fastq.gz | \
      !{params.python} !{params.script_dir}/filter_fastx.py \
      convert fastq fasta > qcs.fasta

    # search for pairs
    !{params.usearch} --usearch_global qcs.fasta \
      --db hqcs.fasta \
      --userout pairfile.txt \
      --userfields query+target \
      --top_hit_only \
      --id !{params.abundance_identity} \
      --maxqt !{params.abundance_max_length_ratio} \
      -qmask none \
      --strand plus \
      --maxaccepts !{params.max_accepts} \
      --maxrejects !{params.max_rejects} \
      --threads !{params.cpus}

    # write abundance file
    !{params.python} !{params.script_dir}/write_abundances.py \
      < pairfile.txt > abundance_file.txt

    # filter out HQCS with 0 abundance
    !{params.python} !{params.script_dir}/filter_fastx.py \
      abundance fasta fasta abundance_file.txt \
      < hqcs.fasta > hqcs.filtered.fasta

    # gzip everything
    rm -f hqcs.fasta
    for i in `find . ! -type l | grep -E "\\.fasta$|\\.fastq$|\\.txt$"`; do gzip "$i" ; done
    '''
}

process merge_timepoints {

    publishDir params.results_dir

    executor 'local'

    input:
    file 'hqcs*.fastq.gz' from hqcs_files.collect()
    file 'abundance*.txt.gz' from abundance_files.collect()

    output:
    file 'all_hqcs.fasta.gz' into merged_hqcs_out

    """
    zcat hqcs*.fastq.gz > merged_hqcs.fasta
    zcat abundance*.txt.gz > merged_abundances.txt

    # add abundances to ids, for evo_history
    ${params.python} ${params.script_dir}/filter_fastx.py \
      add_abundance fasta fasta merged_abundances.txt \
      < merged_hqcs.fasta > all_hqcs.fasta

    # gzip everything
    gzip merged_hqcs.fasta
    gzip all_hqcs.fasta
    gzip merged_abundances.txt
    """
}

/* ************************************************************************** */
/* ALIGNMENT SUB-PIPELINE */

process alignment_pipeline {

    publishDir params.results_dir

    time params.slow_time
    cpus params.cpus
    time params.slow_time

    input:
    file 'hqcs.fasta.gz' from merged_hqcs_out

    output:
    file 'msa.fasta.gz' into msa_out
    file 'msa.aa.fasta.gz' into msa_aa_out

    shell:
    '''
    zcat hqcs.fasta.gz > hqcs.fasta

    !{params.python} !{params.script_dir}/translate.py \
      < hqcs.fasta > hqcs_protein.fasta

    !{params.mafft} --ep 0.5 --quiet --preservecase \
      --thread !{params.cpus} \
      hqcs_protein.fasta > msa.aa.fasta

    !{params.python} !{params.script_dir}/backtranslate.py \
      msa.aa.fasta hqcs.fasta msa.fasta

    # compress everything
    rm -f hqcs.fasta
    for i in `find . ! -type l | grep -E "\\.fasta$|\\.fastq$|\\.txt$"`; do gzip "$i" ; done
    '''
}

/* ************************************************************************** */
/* ANALYSIS SUB-PIPELINE */

process dates_json_task {

    publishDir params.results_dir

    executor 'local'

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

process abundances_json {

    publishDir params.results_dir

    time params.crazy_time

    input:
    file 'msa.fasta.gz' from msa_out

    output:
    file 'abundances.json' into abundances_json_out

    """
    #!${params.python}

    import gzip
    import json
    from Bio import SeqIO
    from flea.util import id_to_abundance

    with gzip.open('msa.fasta.gz', 'rt') as handle:
        records = SeqIO.parse(handle, 'fasta')
        outdict = dict((r.id, id_to_abundance(r.id)) for r in records)
    with open('abundances.json', 'w') as handle:
        json.dump(outdict, handle, separators=(",\\n", ":"))
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
    file 'msa.fasta.gz' from msa_out
    val oldest_label

    output:
    file 'mrca.fasta.gz' into mrca_1, mrca_2, mrca_3
    file 'mrca_translated.fasta.gz' into mrca_translated_1, mrca_translated_2

    shell:
    '''
    zcat msa.fasta.gz | \
    !{params.python} !{params.script_dir}/filter_fastx.py \
      prefix fasta fasta !{oldest_label} \
      > oldest_seqs.fasta

    !{params.python} !{params.script_dir}/DNAcons.py \
      --keep-gaps --codon --id MRCA \
      -o mrca.fasta \
      --abundances \
      oldest_seqs.fasta

    !{params.python} !{params.script_dir}/translate.py --gapped \
      < mrca.fasta > mrca_translated.fasta

    for i in `find . ! -type l | grep -E "\\.fasta$|\\.fastq$|\\.txt$"`; do gzip "$i" ; done
    '''
}

// TODO: why do we have to duplicate outputs here?
process add_mrca {

    executor 'local'

    input:
    file 'mrca.fasta.gz' from mrca_1
    file 'msa.fasta.gz' from msa_out

    output:
    file 'msa_with_mrca.fasta.gz' into msa_with_mrca_1, msa_with_mrca_2

    "zcat mrca.fasta.gz msa.fasta.gz | gzip >  msa_with_mrca.fasta.gz"
}

process infer_tree {

    time params.slow_time

    cpus params.cpus
    time params.slow_time

    input:
    file 'msa.fasta.gz' from msa_with_mrca_1

    output:
    file 'tree.txt' into tree_out

    """
    export OMP_NUM_THREADS=${params.cpus}
    zcat msa.fasta.gz | ${params.fasttree} -gtr -nt > tree.txt
    """
}

process reroot {

    time params.fast_time

    input:
    file 'tree.txt' from tree_out

    output:
    file 'tree.rooted.txt' into rooted_tree_1, rooted_tree_2

    """
    #!${params.python}

    from Bio import Phylo

    tree = next(Phylo.parse('tree.txt', 'newick'))
    clade = next(tree.find_clades('MRCA'))
    tree.root_with_outgroup(clade)

    # also rename for HyPhy
    for i, node in enumerate(tree.get_nonterminals()):
        node.confidence = None
        if node.name != 'MRCA':
            node.name = "ancestor_{}".format(i)
    Phylo.write([tree], 'tree.rooted.txt', 'newick')
    """
}

process tree_json {

    publishDir params.results_dir

    executor 'local'

    input:
    file 'tree.txt' from rooted_tree_1

    output:
    file 'trees.json' into tree_json_out

    """
    #!${params.python}
    import json

    with open('tree.txt') as handle:
        newick_string = handle.read()
    result = {'tree': newick_string}
    with open('trees.json', 'w') as handle:
        json.dump(result, handle, separators=(",\\n", ":"))
    """
}

process js_divergence {

    publishDir params.results_dir

    time params.slow_time

    input:
    file 'msa.aa.fasta.gz' from msa_aa_out
    file 'metadata' from metadata_3

    output:
    file 'js_divergence.json' into js_divergence_json

    """
    zcat msa.aa.fasta.gz > msa.aa.fasta

    ${params.python} ${params.script_dir}/js_divergence.py \
      msa.aa.fasta metadata js_divergence.json

    rm -f msa.aa.fasta
    """
}

process manifold_embedding {

    publishDir params.results_dir

    time params.crazy_time

    input:
    file 'msa.fasta.gz' from msa_out

    output:
    file 'manifold.json' into manifold_json_out

    """
    zcat msa.fasta.gz > msa.fasta
    ${params.usearch} -calc_distmx msa.fasta -distmxout dmatrix.txt \
      -format tabbed_pairs

    ${params.python} ${params.script_dir}/manifold_embed.py \
      --n-jobs 1 dmatrix.txt manifold.json

    rm -f msa.fasta
    gzip dmatrix.txt
    """
}

// TODO: avoid full paths
// TODO: why do command line arguments not work here?
process reconstruct_ancestors {

    time params.slow_time

    input:
    file 'msa.fasta.gz' from msa_with_mrca_2
    file 'msa.aa.fasta.gz' from msa_aa_out
    file 'tree.rooted.txt' from rooted_tree_2

    output:
    file 'msa.aa.ancestors.fasta.gz' into msa_aa_ancestors_out

    shell:
    '''
    zcat msa.fasta.gz > msa.fasta
    zcat msa.aa.fasta.gz > msa.aa.fasta

    echo $(pwd)/msa.fasta >> stdin
    echo $(pwd)/tree.rooted.txt >> stdin
    echo $(pwd)/ancestors.fasta >> stdin
    echo HKY85 >> stdin
    echo 2 >> stdin
    !{params.hyphy} !{params.hyphy_dir}/reconstructAncestors.bf < stdin

    !{params.python} !{params.script_dir}/translate.py --gapped \
      < ancestors.fasta > ancestors.aa.fasta

    cat msa.aa.fasta ancestors.aa.fasta > 'msa.aa.ancestors.fasta'

    # compress all
    rm -f msa.fasta msa.aa.fasta
    for i in `find . ! -type l | grep -E "\\.fasta$|\\.fastq$|\\.txt$"`; do gzip "$i" ; done
    '''
}

process coordinates_json {

    publishDir params.results_dir

    time params.fast_time

    cpus params.cpus
    time params.slow_time

    input:
    file 'mrca.aa.fasta.gz' from mrca_translated_1

    output:
    file 'coordinates.json' into coordinates_json_out

    """
    zcat mrca.aa.fasta.gz > mrca.aa.fasta
    cat mrca.aa.fasta ${params.reference_protein} > pair.fasta

    ${params.mafft} --ep 0.5 --quiet --preservecase \
      --thread ${params.cpus} \
      pair.fasta > aligned.fasta

    ${params.python} ${params.script_dir}/coordinates_json.py \
      mrca.aa.fasta aligned.fasta ${params.reference_coordinates} coordinates.json

    rm -f mrca.aa.fasta
    """
}

process sequences_json {
    publishDir params.results_dir

    time params.fast_time

    input:
    file 'msa.fasta.gz' from msa_aa_ancestors_out
    file 'mrca.fasta.gz' from mrca_translated_2
    file 'coordinates.json' from coordinates_json_out
    file 'metadata' from metadata_4

    output:
    file 'sequences.json' into sequences_json_out

    """
    zcat msa.fasta.gz > msa.fasta
    zcat mrca.fasta.gz > mrca.fasta
    ${params.python} ${params.script_dir}/sequences_json.py \
      msa.fasta mrca.fasta coordinates.json metadata \
      ${params.reference_protein} ${params.reference_coordinates} \
      sequences.json

    rm -f msa.fasta mrca.fasta
    """
}

process replace_stop_codons {

    time params.fast_time

    input:
    file 'msa.fasta.gz' from msa_out

    output:
    file 'msa.no_stops.fasta.gz' into msa_no_stops

    """
    zcat msa.fasta.gz |
    ${params.python} ${params.script_dir}/filter_fastx.py \
      stop_codons fasta fasta | gzip > msa.no_stops.fasta.gz
    """
}

// TODO: why do we need to split output here, but not elsewhere?
process seq_dates {

    time params.fast_time

    input:
    file 'msa.fasta.gz' from msa_out
    file 'metadata' from metadata_5

    output:
    file 'dates.json' into seq_dates_1, seq_dates_2

    """
    #!${params.python}

    import gzip
    import json
    from Bio import SeqIO
    from flea.util import get_date_dict
    from flea.util import id_to_label

    date_dict = get_date_dict('metadata')
    with gzip.open('msa.fasta.gz', 'rt') as handle:
        records = SeqIO.parse(handle, "fasta")
        outdict = dict((r.id, date_dict[id_to_label(r.id)]) for r in records)
    with open('dates.json', "w") as handle:
        json.dump(outdict, handle, separators=(",\\n", ":"))
    """
}

process region_coords {

    publishDir params.results_dir

    time params.fast_time

    input:
    file 'mrca.fasta.gz' from mrca_2

    output:
    file 'region_coords.json' into region_coords_json

    shell:
    '''
    zcat mrca.fasta.gz > mrca.fasta
    !{params.hyphy} !{params.hyphy_dir}/HXB2partsSplitter.bf \
      $(pwd)/mrca.fasta $(pwd)/region_coords.json
    rm -f mrca.fasta
    '''
}


process evo_history {

    time params.crazy_time

    input:
    file 'msa.no_stops.fasta.gz' from msa_no_stops
    file 'dates.json' from seq_dates_1
    file 'region_coords.json' from region_coords_json

    output:
    file 'rates_pheno.tsv' into rates_pheno

    when:
    params.do_evo_history

    shell:
    '''
    zcat msa.no_stops.fasta.gz > msa.no_stops.fasta

    !{params.hyphy} !{params.hyphy_dir}/obtainEvolutionaryHistory.bf \
      $(pwd)/msa.no_stops.fasta $(pwd)/dates.json $(pwd)/region_coords.json $(pwd)/rates_pheno.tsv

    rm -f msa.no_stops.fasta
    '''
}

process rates_pheno_json {

    publishDir params.results_dir

    time params.fast_time

    input:
    file 'rates_pheno.tsv' from rates_pheno

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

    publishDir params.results_dir

    cpus params.crazy_cpus
    time params.crazy_time

    input:
    file 'msa.no_stops.fasta.gz' from msa_no_stops
    file 'dates.json' from seq_dates_2
    file 'mrca.fasta.gz' from mrca_3

    output:
    file 'rates.json' into rates_json

    when:
    params.do_fubar

    shell:
    '''
    zcat msa.no_stops.fasta.gz > msa.no_stops.fasta
    zcat mrca.fasta.gz > mrca.fasta

    !{params.hyphympi} !{params.hyphy_dir}/runFUBAR.bf \
      $(pwd)/msa.no_stops.fasta $(pwd)/dates.json $(pwd)/mrca.fasta $(pwd)/ $(pwd)/rates.json

    rm -f msa.no_stops.fasta mrca.fasta
    '''
}

/* ************************************************************************** */
/* DIAGNOSIS SUB-PIPELINE */

process diagnose {

    publishDir params.results_dir

    cpus params.cpus
    time params.crazy_time

    input:
    file 'qcs*.fastq.gz' from qcs_final_4.map{ it[0] }.collect()
    file 'hqcs.fasta.gz' from merged_hqcs_out
    file 'hqcs.msa.fasta.gz' from msa_out
    file 'hqcs.msa.aa.fasta.gz' from msa_aa_out

    output:
    file 'diagnosis_results/*' into diagnosis_results

    when:
    params.do_diagnosis

    shell:
    '''
    zcat hqcs.fasta.gz > hqcs.fasta
    zcat hqcs.msa.fasta.gz > hqcs.msa.fasta
    zcat hqcs.msa.aa.fasta.gz > hqcs.msa.aa.fasta

    # convert to fasta
    zcat qcs*.fastq.gz |
      !{params.python} !{params.script_dir}/filter_fastx.py \
        convert fastq fasta > qcs.fasta

    # search for pairs
    !{params.usearch} --usearch_global qcs.fasta \
      --db hqcs.fasta \
      --userout pairfile.txt \
      --userfields query+target \
      --top_hit_only \
      --id !{params.qcs_to_hqcs_identity} \
      --maxqt !{params.max_qcs_length_ratio} \
      --qmask none \
      --strand plus \
      --maxaccepts !{params.max_accepts} \
      --maxrejects !{params.max_rejects} \
      --threads !{params.cpus}

    # combine pairs
    mkdir alignments
    !{params.python} !{params.script_dir}/pairs_to_fasta.py \
      qcs.fasta hqcs.fasta pairfile.txt alignments

    # align and insert gaps
    !{params.parallel} -j !{params.cpus} \
      '!{params.bealign} -a codon {} {}.bam' ::: alignments/*unaligned.fasta

    !{params.parallel} -j !{params.cpus} \
      '!{params.bam2msa} {} {}.fasta' ::: alignments/*.bam

    !{params.parallel} -j !{params.cpus} \
      '!{params.python} !{params.script_dir}/insert_gaps.py {} hqcs.msa.fasta {}.gapped' ::: alignments/*.bam.fasta

    # merge all
    cat alignments/*gapped > qcs.msa.fasta

    # replace gapped codons
    !{params.python} !{params.script_dir}/filter_fastx.py \
      gap_codons fasta fasta < qcs.msa.fasta > qcs.msa.nogaps.fasta

    !{params.python} !{params.script_dir}/translate.py --gapped \
      < qcs.msa.nogaps.fasta > qcs.msa.aa.fasta

    # run diagnosis
    mkdir diagnosis_results
    !{params.python} !{params.script_dir}/diagnose.py \
      hqcs.msa.aa.fasta qcs.msa.aa.fasta diagnosis_results

    # gzip all
    rm -f hqcs.msa.fasta hqcs.msa.aa.fasta hqcs.fasta
    for i in `find . ! -type l | grep -E "\\.fasta$|\\.fastq$|\\.txt$"`; do gzip "$i" ; done
    '''
}