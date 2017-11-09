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


// TODO: how launch with web ui and monitor progress
// TODO: progressively retry with longer times, if there is a timeout
// TODO: update tests
// TODO: tune maxaccepts and maxrejects
// TODO: combine all time points for inframe db for frame correction


params.infile = "$HOME/flea/data/P018/data/metadata"
params.results_dir = "results"


// TODO: how to avoid duplicating?
Channel.fromPath(params.infile)
    .into { metadata_1; metadata_2; metadata_3; metadata_4; metadata_5; metadata_6 }


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


compress_cmd = """for i in `find . ! -type l | grep -E "\\.fasta\$|\\.fastq\$|\\.txt\$|\\.dst\$"`; do gzip "\$i" ; done"""


// TODO: train head/tail HMM on all sequences from all time points
hmm_train_flag = (params.train_hmm ? '--train' : '')
process quality_pipeline {

    tag { label }

    publishDir params.results_dir

    time params.slow_time

    afterScript compress_cmd

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
    !{params.python} !{workflow.projectDir}/flea/trim_tails.py \
      --n-jobs !{params.cpus} --fastq \
      !{hmm_train_flag} --max-iters !{params.train_hmm_max_iters} \
      qfiltered.fastq trimmed.fastq

    # filter runs
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      runs fastq fasta !{params.run_length} < trimmed.fastq > no_runs.fasta

    # filter against contaminants
    !{params.usearch} --usearch_global no_runs.fasta \
      --db !{params.contaminants_db} \
      --notmatched uncontam.fasta \
      --matched contam.fasta \
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
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      userout fastq fastq userout.txt \
      < trimmed.fastq > filtered.fastq

    # length filter
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      length fastq fastq \
      !{minlen} !{maxlen} \
      < filtered.fastq > !{label}.qcs.fastq
    '''
}

/* ************************************************************************** */
/* CONSENSUS SUB-PIPELINE */


process cluster {

    tag { label }

    publishDir params.results_dir

    time params.slow_time

    input:
    set 'qcs.fastq.gz', label from qcs_final_1

    output:
    set '*.clusters.uc', label into cluster_out

    shell:
    '''
    zcat qcs.fastq.gz > qcs.fastq

    # sort by error rate
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      !{params.before_cluster} fastq fastq \
      < qcs.fastq > qcs.sorted.fastq

    # cluster
    !{params.usearch} --cluster_fast qcs.sorted.fastq \
      -uc !{label}.clusters.uc \
      --id !{params.cluster_identity} \
      --minsl !{params.min_length_ratio} \
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

    time params.slow_time

    afterScript compress_cmd

    input:
    set 'clusters.uc', 'qcs.fastq.gz', label from consensus_input

    output:
    set '*.clusters.consensus.fasta.gz', label into consensus_out

    shell:
    '''
    zcat qcs.fastq.gz | \
    !{params.python} !{workflow.projectDir}/flea/cluster_fastq.py \
      --minsize !{params.min_cluster_size} \
      clusters.uc .

    # function sample clusters and do mafft consensus
    doconsensus() {
        !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
          sample fastq fasta \
          !{params.min_cluster_size} !{params.max_cluster_size} \
          < ${1} > ${1}.sampled.fasta

        !{params.mafft} --ep 0.5 --quiet --preservecase \
          ${1}.sampled.fasta > ${1}.sampled.aligned.fasta

        # get cluster number so we can put it in the record id
        number=$(echo ${1} | cut -d '_' -f 2)

        !{params.python} !{workflow.projectDir}/flea/DNAcons.py \
          --seed 0 \
          -o ${1}.consensus.fasta \
          --name !{label}_consensus_${number} \
          ${1}.sampled.aligned.fasta
    }
    export -f doconsensus

    # run in parallel
    !{params.parallel} -j !{params.cpus} 'doconsensus {}' ::: *_raw.fastq

    cat *.consensus.fasta > !{label}.clusters.consensus.fasta

    # check that all sequences are present
    n_expected=`ls *_raw.fastq | wc -l`
    n_found=`grep ">" !{label}.clusters.consensus.fasta | wc -l`
    if [ "$n_expected" -ne "$n_found" ]; then
        echo "ERROR: some consensus sequences are missing"
        exit 1
    fi

    rm -f qcs.fastq
    '''
}


allow_stop_codons = params.do_frame_correction ? "false" : "true"

process inframe_unique_hqcs {

    tag { label }

    publishDir params.results_dir

    afterScript compress_cmd

    input:
    set 'consensus.fasta.gz', label from consensus_out

    output:
    set '*.consensus.fasta.gz', '*.consensus.unique.fasta.gz', label into inframe_unique_out_1,
        inframe_unique_out_2,
        inframe_unique_out_3

    shell:
    '''
    zcat consensus.fasta.gz > consensus.fasta

    # inframe
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      inframe fasta fasta !{allow_stop_codons} \
      < consensus.fasta > consensus.inframe.fasta

    # unique
    !{params.usearch} --fastx_uniques consensus.inframe.fasta \
      --fastaout !{label}.consensus.unique.fasta \
      --threads !{params.cpus}

    cp consensus.fasta !{label}.consensus.fasta

    rm -f consensus.fasta
    '''
}


inframe_unique_out_1
  .map { it -> it[1] }
  .set { inframe_unique_dbs }


process make_inframe_db {

    publishDir params.results_dir

    afterScript compress_cmd

    when:
    params.do_frame_correction

    input:
    file '*.consensus.unique.fasta.gz' from inframe_unique_dbs.collect()

    output:
    file 'inframedb.fasta.gz' into inframe_db_out

    shell:
    '''
    zcat *.consensus.unique.fasta.gz > inframedb.fasta
    '''
}

process frame_correction {

    tag { label }

    publishDir params.results_dir

    time params.slow_time

    afterScript compress_cmd

    when:
    params.do_frame_correction

    input:
    set 'consensus.fasta.gz', 'unused.db.fasta.gz', label from inframe_unique_out_2
    file 'inframedb.fasta.gz' from inframe_db_out

    output:
    set '*.consensus.unique.corrected.fasta.gz', label into frame_correction_out

    shell:
    '''
    zcat consensus.fasta.gz > consensus.fasta
    zcat inframedb.fasta.gz > inframedb.fasta

    # search
    !{params.usearch} --usearch_global consensus.fasta \
      --db inframedb.fasta \
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

    # frame correction
    !{params.python} !{workflow.projectDir}/flea/frame_correction.py \
      --deletion-strategy=reference \
      --calns=calnfile.txt \
      pairfile.fasta corrected.fasta

    # filter inframe
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      inframe fasta fasta true < corrected.fasta > corrected.inframe.fasta

    # deduplicate
    !{params.usearch} --fastx_uniques corrected.inframe.fasta \
      --fastaout !{label}.consensus.unique.corrected.fasta \
      --threads !{params.cpus}

    rm -f consensus.fasta
    rm -f consensus.db.fasta
    '''
}


// NOTE: if previous steps are cached for both values of
// `do_frame_correction`, changing the value of `do_frame_correction`
// doesn't update the symlinks downstream.
if( params.do_frame_correction ) {
    frame_correction_out
      .phase (qcs_final_3) { it[1] }
      .map { [ it[0][0], it[1][0], it[0][1] ] }
      .set { compute_copynumbers_input }
} else {
    inframe_unique_out_3
      .phase (qcs_final_3) { it[-1] }
      .map { [ it[0][1], it[1][0], it[0][2] ] }
      .set { compute_copynumbers_input }
}


process compute_copynumbers {

    tag { label }

    time params.slow_time

    afterScript compress_cmd

    input:
    set 'hqcs.fasta.gz', 'qcs.fastq.gz', label from compute_copynumbers_input

    output:
    file 'hqcs.filtered.fasta.gz' into hqcs_files
    file 'copynumber_file.txt.gz' into copynumber_files

    shell:
    '''
    zcat hqcs.fasta.gz > hqcs.fasta

    # convert to fasta for usearch
    zcat qcs.fastq.gz | \
      !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      convert fastq fasta > qcs.fasta

    # search for pairs
    !{params.usearch} --usearch_global qcs.fasta \
      --db hqcs.fasta \
      --userout pairfile.txt \
      --userfields query+target \
      --top_hit_only \
      --id !{params.copynumber_identity} \
      --maxqt !{params.copynumber_max_length_ratio} \
      -qmask none \
      --strand plus \
      --maxaccepts !{params.max_accepts} \
      --maxrejects !{params.max_rejects} \
      --threads !{params.cpus}

    # write copynumber file
    !{params.python} !{workflow.projectDir}/flea/write_copynumbers.py \
      < pairfile.txt > copynumber_file.txt

    # filter out HQCS with 0 copynumber
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      copynumber fasta fasta copynumber_file.txt \
      < hqcs.fasta > hqcs.filtered.fasta

    rm -f hqcs.fasta
    '''
}

process merge_timepoints {

    publishDir params.results_dir

    afterScript compress_cmd

    executor 'local'
    cpus 1

    input:
    file 'hqcs*.fastq.gz' from hqcs_files.collect()
    file 'copynumber*.txt.gz' from copynumber_files.collect()

    output:
    file 'all_hqcs.fasta.gz' into merged_hqcs_out

    """
    zcat hqcs*.fastq.gz > merged_hqcs.fasta
    zcat copynumber*.txt.gz > merged_copynumbers.txt

    # add copynumbers to ids, for evo_history
    ${params.python} ${workflow.projectDir}/flea/filter_fastx.py \
      add_copynumber fasta fasta merged_copynumbers.txt \
      < merged_hqcs.fasta > all_hqcs.fasta
    """
}

/* ************************************************************************** */
/* ALIGNMENT SUB-PIPELINE */

process alignment_pipeline {

    publishDir params.results_dir

    time params.slow_time

    afterScript compress_cmd

    input:
    file 'hqcs.fasta.gz' from merged_hqcs_out

    output:
    file 'msa.fasta.gz' into msa_out, msa_out_2
    file 'msa.aa.fasta.gz' into msa_aa_out

    shell:
    '''
    zcat hqcs.fasta.gz > hqcs.fasta

    !{params.python} !{workflow.projectDir}/flea/translate.py \
      < hqcs.fasta > hqcs_protein.fasta

    !{params.mafft} --ep 0.5 --quiet --preservecase \
      --thread !{params.cpus} \
      hqcs_protein.fasta > msa.aa.fasta

    !{params.python} !{workflow.projectDir}/flea/backtranslate.py \
      msa.aa.fasta hqcs.fasta msa.fasta

    rm -f hqcs.fasta
    '''
}

/* ************************************************************************** */
/* ANALYSIS SUB-PIPELINE */

process dates_json_task {

    publishDir params.results_dir

    executor 'local'
    cpus 1

    when:
    params.do_analysis

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

    publishDir params.results_dir

    executor 'local'
    cpus 1

    when:
    params.do_analysis

    input:
    file 'msa.fasta.gz' from msa_out

    output:
    file 'copynumbers.json' into copynumbers_json_out

    """
    #!${params.python}

    import gzip
    import json
    from Bio import SeqIO
    from flea.util import id_to_copynumber

    with gzip.open('msa.fasta.gz', 'rt') as handle:
        records = SeqIO.parse(handle, 'fasta')
        outdict = dict((r.id, id_to_copynumber(r.id)) for r in records)
    with open('copynumbers.json', 'w') as handle:
        json.dump(outdict, handle, separators=(",\\n", ":"))
    """
}

process get_oldest_label {

    executor 'local'
    cpus 1

    when:
    params.do_analysis

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

    afterScript compress_cmd

    when:
    params.do_analysis

    input:
    file 'msa.fasta.gz' from msa_out
    val oldest_label

    output:
    file 'mrca.fasta.gz' into mrca_1, mrca_2, mrca_3, mrca_4
    file 'mrca_translated.fasta.gz' into mrca_translated_1, mrca_translated_2, mrca_translated_3

    shell:
    '''
    zcat msa.fasta.gz | \
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      prefix fasta fasta !{oldest_label} \
      > oldest_seqs.fasta

    !{params.python} !{workflow.projectDir}/flea/DNAcons.py \
      --seed 0 \
      --keep-gaps \
      --codon \
      --copynumbers \
      --name MRCA \
      -o mrca.fasta \
      oldest_seqs.fasta

    !{params.python} !{workflow.projectDir}/flea/translate.py --gapped \
      < mrca.fasta > mrca_translated.fasta
    '''
}

// TODO: why do we have to duplicate outputs here?
process add_mrca {

    executor 'local'
    cpus 1

    when:
    params.do_analysis

    input:
    file 'mrca.fasta.gz' from mrca_1
    file 'msa.fasta.gz' from msa_out

    output:
    file 'msa_with_mrca.fasta.gz' into msa_with_mrca_1, msa_with_mrca_2

    "zcat mrca.fasta.gz msa.fasta.gz | gzip >  msa_with_mrca.fasta.gz"
}

process infer_tree {

    time params.slow_time

    when:
    params.do_analysis

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

    executor 'local'
    cpus 1

    when:
    params.do_analysis

    input:
    file 'tree.txt' from tree_out

    output:
    file 'tree.rooted.txt' into rooted_tree_1, rooted_tree_2, rooted_tree_3

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
    cpus 1

    when:
    params.do_analysis

    input:
    file 'tree.txt' from rooted_tree_1

    output:
    file 'trees.json' into trees_json_out

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

    when:
    params.do_analysis

    input:
    file 'msa.aa.fasta.gz' from msa_aa_out
    file 'mrca.aa.fasta.gz' from mrca_translated_1
    file 'metadata' from metadata_3

    output:
    file 'js_divergence.json' into js_divergence_json_out

    """
    zcat msa.aa.fasta.gz > msa.aa.fasta
    zcat mrca.aa.fasta.gz > mrca.aa.fasta

    ${params.python} ${workflow.projectDir}/flea/js_divergence.py \
      msa.aa.fasta mrca.aa.fasta metadata js_divergence.json

    rm -f msa.aa.fasta mrca.aa.fasta
    """
}

process manifold_embedding {

    publishDir params.results_dir

    time params.slow_time

    afterScript compress_cmd

    when:
    params.do_analysis

    input:
    file 'msa.fasta.gz' from msa_out
    file 'metadata' from metadata_4

    output:
    file 'manifold.json' into manifold_json_out

    """
    zcat msa.fasta.gz > msa.fasta
    ${params.tn93} -t ${params.tn93_threshold} -o dmatrix.dst msa.fasta

    ${params.python} ${workflow.projectDir}/flea/manifold_embed.py \
      --n-jobs 1 dmatrix.dst metadata manifold.json

    rm -f msa.fasta
    """
}

// TODO: avoid full paths
// TODO: why do command line arguments not work here?
process reconstruct_ancestors {

    time params.slow_time

    afterScript compress_cmd

    when:
    params.do_analysis

    input:
    file 'msa.fasta.gz' from msa_with_mrca_2
    file 'msa.aa.fasta.gz' from msa_aa_out
    file 'tree.rooted.txt' from rooted_tree_2

    output:
    file 'msa.aa.ancestors.fasta.gz' into msa_aa_ancestors_out
    file 'ancestors.fasta.gz' into ancestors_out

    shell:
    '''
    zcat msa.fasta.gz > msa.fasta
    zcat msa.aa.fasta.gz > msa.aa.fasta

    echo $(pwd)/msa.fasta >> stdin
    echo $(pwd)/tree.rooted.txt >> stdin
    echo $(pwd)/ancestors.fasta >> stdin
    echo HKY85 >> stdin
    echo 2 >> stdin
    !{params.hyphy} !{workflow.projectDir}/hyphy_scripts/reconstructAncestors.bf < stdin

    !{params.python} !{workflow.projectDir}/flea/translate.py --gapped \
      < ancestors.fasta > ancestors.aa.fasta

    cat msa.aa.fasta ancestors.aa.fasta > 'msa.aa.ancestors.fasta'

    rm -f msa.fasta msa.aa.fasta
    '''
}

process coordinates_json {

    publishDir params.results_dir

    when:
    params.do_analysis

    input:
    file 'mrca.aa.fasta.gz' from mrca_translated_2

    output:
    file 'coordinates.json' into coordinates_json_out_1, coordinates_json_out_2

    """
    zcat mrca.aa.fasta.gz > mrca.aa.fasta
    cat mrca.aa.fasta ${params.reference_protein} > pair.fasta

    ${params.mafft} --ep 0.5 --quiet --preservecase \
      --thread ${params.cpus} \
      pair.fasta > aligned.fasta

    ${params.python} ${workflow.projectDir}/flea/coordinates_json.py \
      mrca.aa.fasta aligned.fasta ${params.reference_coordinates} coordinates.json

    rm -f mrca.aa.fasta
    """
}

process sequences_json {
    publishDir params.results_dir

    when:
    params.do_analysis

    input:
    file 'msa.fasta.gz' from msa_aa_ancestors_out
    file 'mrca.fasta.gz' from mrca_translated_3
    file 'coordinates.json' from coordinates_json_out_1
    file 'metadata' from metadata_5

    output:
    file 'sequences.json' into sequences_json_out

    """
    zcat msa.fasta.gz > msa.fasta
    zcat mrca.fasta.gz > mrca.fasta
    ${params.python} ${workflow.projectDir}/flea/sequences_json.py \
      msa.fasta mrca.fasta coordinates.json metadata \
      ${params.reference_protein} ${params.reference_coordinates} \
      sequences.json

    rm -f msa.fasta mrca.fasta
    """
}

process replace_stop_codons {

    when:
    params.do_analysis

    input:
    file 'msa.fasta.gz' from msa_out

    output:
    file 'msa.no_stops.fasta.gz' into msa_no_stops

    """
    zcat msa.fasta.gz |
    ${params.python} ${workflow.projectDir}/flea/filter_fastx.py \
      stop_codons fasta fasta | gzip > msa.no_stops.fasta.gz
    """
}

// TODO: why do we need to split output here, but not elsewhere?
process seq_dates {

    when:
    params.do_analysis

    input:
    file 'msa.fasta.gz' from msa_out
    file 'metadata' from metadata_6

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

    when:
    params.do_analysis

    input:
    file 'mrca.fasta.gz' from mrca_2

    output:
    file 'region_coords.json' into region_coords_json

    shell:
    '''
    zcat mrca.fasta.gz > mrca.fasta
    !{params.hyphy} !{workflow.projectDir}/hyphy_scripts/HXB2partsSplitter.bf \
      $(pwd)/mrca.fasta $(pwd)/region_coords.json
    rm -f mrca.fasta
    '''
}


process evo_history {

    publishDir params.results_dir

    // need to make module command available on compute node
    beforeScript params.module_script
    module params.modules

    time params.slow_time

    when:
    params.do_analysis && params.do_evo_history

    input:
    file 'msa.no_stops.fasta.gz' from msa_no_stops
    file 'dates.json' from seq_dates_1
    file 'region_coords.json' from region_coords_json

    output:
    file 'rates_pheno.tsv' into rates_pheno

    shell:
    '''
    zcat msa.no_stops.fasta.gz > msa.no_stops.fasta

    !{params.hyphympi} !{workflow.projectDir}/hyphy_scripts/obtainEvolutionaryHistory.bf \
      $(pwd)/msa.no_stops.fasta $(pwd)/dates.json $(pwd)/region_coords.json $(pwd)/rates_pheno.tsv

    rm -f msa.no_stops.fasta
    '''
}

process rates_pheno_json {

    publishDir params.results_dir

    when:
    params.do_analysis && params.do_evo_history

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

    // need to make module command available on compute node
    beforeScript params.module_script
    module params.modules

    time params.crazy_time

    when:
    params.do_analysis && params.do_fubar

    input:
    file 'msa.no_stops.fasta.gz' from msa_no_stops
    file 'dates.json' from seq_dates_2
    file 'mrca.fasta.gz' from mrca_3

    output:
    file 'rates.json' into rates_json_out

    shell:
    '''
    zcat msa.no_stops.fasta.gz > msa.no_stops.fasta
    zcat mrca.fasta.gz > mrca.fasta

    !{params.hyphympi} !{workflow.projectDir}/hyphy_scripts/runFUBAR.bf \
      $(pwd)/msa.no_stops.fasta $(pwd)/dates.json $(pwd)/mrca.fasta $(pwd)/ $(pwd)/rates.json

    rm -f msa.no_stops.fasta mrca.fasta
    '''
}

Channel.empty()
        .mix(
             coordinates_json_out_2,
             copynumbers_json_out,
             dates_json_out,
             js_divergence_json_out,
             manifold_json_out,
             sequences_json_out,
             trees_json_out,
             rates_pheno_json_out,
             rates_json_out
             )
  .flatten()
  .collect()
  .set { json_inputs }

process combine_results {
    publishDir params.results_dir

    executor 'local'
    cpus 1

    input:
    file '*' from json_inputs

    output:
    file 'session.json' into session_json

    """
    ${params.python} ${workflow.projectDir}/flea/convert_jsons.py .
    """
}

Channel.empty()
        .mix(
             ancestors_out,
	     mrca_4,
             msa_out_2,
             rooted_tree_3,
             )
  .flatten()
  .collect()
  .set { zip_inputs }

process zip_results {
    publishDir params.results_dir

    executor 'local'
    cpus 1

    input:
    file '*' from zip_inputs

    output:
    file 'session.zip' into session_zip

    """
    zcat msa.fasta.gz > msa.fasta
    zcat mrca.fasta.gz > mrca.fasta
    zcat ancestors.fasta.gz > ancestors.fasta

    zip session.zip *.txt *.fasta

    rm -rf *.fasta
    """
}



/* ************************************************************************** */
/* DIAGNOSIS SUB-PIPELINE */

process diagnose {

    publishDir params.results_dir

    time params.crazy_time

    afterScript compress_cmd

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
      !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
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
    !{params.python} !{workflow.projectDir}/flea/pairs_to_fasta.py \
      qcs.fasta hqcs.fasta pairfile.txt alignments

    # align and insert gaps
    !{params.parallel} -j !{params.cpus} \
      '!{params.bealign} -a codon {} {}.bam' ::: alignments/*unaligned.fasta

    !{params.parallel} -j !{params.cpus} \
      '!{params.bam2msa} {} {}.fasta' ::: alignments/*.bam

    !{params.parallel} -j !{params.cpus} \
      '!{params.python} !{workflow.projectDir}/flea/insert_gaps.py {} hqcs.msa.fasta {}.gapped' ::: alignments/*.bam.fasta

    # merge all
    cat alignments/*gapped > qcs.msa.fasta

    # replace gapped codons
    !{params.python} !{workflow.projectDir}/flea/filter_fastx.py \
      gap_codons fasta fasta < qcs.msa.fasta > qcs.msa.nogaps.fasta

    !{params.python} !{workflow.projectDir}/flea/translate.py --gapped \
      < qcs.msa.nogaps.fasta > qcs.msa.aa.fasta

    # run diagnosis
    mkdir diagnosis_results
    !{params.python} !{workflow.projectDir}/flea/diagnose.py \
      hqcs.msa.aa.fasta qcs.msa.aa.fasta diagnosis_results

    rm -f hqcs.msa.fasta hqcs.msa.aa.fasta hqcs.fasta
    '''
}
