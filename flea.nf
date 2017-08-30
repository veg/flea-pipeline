#!/usr/bin/env nextflow

/*
 * The FLEA pipeline
 *
 */

// TODO: run on TORQUE
// TODO: set cpus and threads for each task
// TODO: write consensus, alignment, analysis, diagnosis pipelines
// TODO: allow starting from aligned inputs
// TODO: do not hardcode min/max lengths
// TODO: sub-pipelines
// TODO: nextflow repo format
// TODO: infer singleton proc

params.infile = "$HOME/flea/data/P018/data/metadata"

// TODO: how to avoid duplicating?
Channel.fromPath(params.infile)
    .into { metadata_1; metadata_2; metadata_3 }


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
    set 'lenfiltered', label into qcs_final_1, qcs_final_2, qcs_final_3

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
    file 'merged_hqcs' into merged_hqcs_out

    """
    cat hqcs* > merged_hqcs
    """
}

process merge_copynumbers {
    input:
    file 'cn*' from cnfiles.collect()

    output:
    file all_cns into all_cns_out

    """
    cat cn* > all_cns
    """
}

/* ************************************************************************** */
/* ALIGNMENT SUB-PIPELINE */

process translate_hqcs {
    input:
    file 'hqcs' from merged_hqcs_out

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
    file 'dna' from merged_hqcs_out
    file 'aa' from hqcs_aligned

    output:
    file 'msa' into msa_out

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

process dates_json_task {
    input:
    file 'metadata' from metadata_2

    output:
    file 'dates.json' into dates_json_out

    """
    #!${params.python}
    import json
    from flea.util import get_date_dict

    result = get_date_dict('metadata')
    with open('dates.json', 'w') as handle:
        json.dump(result, handle, separators=(",\\n", ":"))
    """
}

// processs copynumber_json {
//     input:
//     file 'cns' from all_cns_out

//     output:
//     file 'copynumbers.json' into copynumbers_json_out

//     """
//     #!${params.python}

//     from flea.utils import id_with_cn

//     d = parse_copynumbers('cns')
//     # add copynumber to name, to match rest of this pipeline
//     result = dict((id_with_cn(key, value), value)
//                   for key, value in d.items())
//     with open('copynumbers.json', 'w') as handle:
//         json.dump(result, handle, separators=(",\\n", ":"))
//     """
// }

// TODO: rewrite as filter_fastx with id prefix
process oldest_seqs {
    input:
    file 'msa' from msa_out

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
    file 'cns' from all_cns_out

    output:
    file 'mrca' into mrca_out

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
    file 'msa' from msa_out
    file 'cns' from all_cns_out

    output:
    file 'msa_with_cn' into msa_id_out

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
    file 'mrca' from mrca_out
    file 'msa' from msa_id_out

    output:
    file 'msa_with_mrca' into msa_with_mrca_out

    "cat mrca msa > msa_with_mrca"
}

process fasttree {
    input:
    file 'msa' from msa_with_mrca_out

    output:
    file tree

    "${params.fasttree} -gtr -nt msa > tree"
}

process reroot {
    input:
    file tree

    output:
    file 'rooted_tree' into rooted_tree_out

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
    input:
    file 'tree' from rooted_tree_out

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

process translate_msa_with_cn {
    input:
    file 'msa' from msa_id_out

    output:
    file msa_translated

    """
    ${params.python} ${params.script_dir}/translate.py --gapped \
      < msa > msa_translated
    """
}

process translate_mrca {
    input:
    file 'mrca' from mrca_out

    output:
    file mrca_translated

    """
    ${params.python} ${params.script_dir}/translate.py --gapped \
      < mrca > mrca_translated
    """
}

process js_divergence {
    input:
    file 'msa' from msa_id_out
    file 'metadata' from metadata_1

    output:
    file 'js_divergence.json' into js_divergence_json

    """
    ${params.python} ${params.script_dir}/js_divergence.py \
      msa metadata js_divergence.json
    """
}

process distance_matrix {
    input:
    file 'msa' from msa_id_out

    output:
    file dmatrix

    """
    ${params.usearch} -calc_distmx msa -distmxout dmatrix \
      -format tabbed_pairs
    """
}

// TODO: threading
process manifold_embedding {
    input:
    file dmatrix

    output:
    file 'manifold.json' into manifold_json_out

    """
    ${params.python} ${params.script_dir}/manifold_embed.py \
      --n-jobs 1 --flip \
      dmatrix manifold.json
    """
}

// TODO: avoid full paths
// TODO: why do command line arguments not work here?
process reconstruct_ancestors {
    input:
    file 'msa' from msa_with_mrca_out
    file 'tree' from rooted_tree_out

    output:
    file 'ancestors' into ancestors_out

    shell:
    '''
    echo $(pwd)/msa >> stdin
    echo $(pwd)/tree >> stdin
    echo $(pwd)/ancestors >> stdin
    echo HKY85 >> stdin
    echo 2 >> stdin
    !{params.hyphy} !{params.hyphy_dir}/reconstructAncestors.bf < stdin
    '''
}

process translate_ancestors {
    input:
    file 'ancestors' from ancestors_out

    output:
    file 'translated' into ancestors_translated_out

    """
    ${params.python} ${params.script_dir}/translate.py --gapped \
      < ancestors > translated
    """
}

process replace_stop_codons {
    input:
    file 'msa' from msa_id_out

    output:
    file 'no_stops' into msa_no_stops

    """
    ${params.python} ${params.script_dir}/filter_fastx.py \
      stop_codons fasta fasta < msa > no_stops
    """
}

// TODO: why do we need to split output here, but not elsewhere?
process seq_dates {
    input:
    file 'msa' from msa_id_out
    file 'metadata' from metadata_3

    output:
    file 'dates' into seq_dates_1, seq_dates_2

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

// FIXME: why does this segfault???

/*
process region_coords {
    input:
    file 'mrca' from mrca_out

    output:
    file 'region_coords.json' into region_coords_json

    shell:
    '''
    !{params.hyphy} !{params.hyphy_dir}/HXB2partsSplitter.bf \
      $(pwd)/mrca $(pwd)/region_coords.json
    '''
}

process evo_history {
    input:
    file 'no_stops' from msa_no_stops
    file 'dates' from seq_dates_1
    file 'region_coords' from region_coords_json

    output:
    file rates_pheno

    shell:
    '''
    !{params.hyphy} !{params.hyphy_dir}/obtainEvolutionaryHistory.bf \
      $(pwd)/no_stops $(pwd)/dates $(pwd)/region_coords
    '''
}

process rates_pheno_json {
    input:
    file rates_pheno

    output:
    file 'rates_pheno.json' into rates_pheno_json_out

    """
    #!${params.python}

    import json

    with open('rates_pheno') as h:
        lines = list(csv.reader(h, delimiter='\t'))
    result = {}
    keys, rest = lines[0], lines[1:]
    result = list(dict((key, value) for key, value in zip(keys, line))
                  for line in rest)
    with open('rates_pheno.json', 'w') as h:
        json.dump(result, h, separators=(",\\n", ":"))
    """
}
*/

process fubar {
    input:
    file 'no_stops' from msa_no_stops
    file 'dates' from seq_dates_2
    file 'mrca' from mrca_out

    output:
    file 'rates.json' into rates_json

    shell:
    '''
    !{params.hyphy} !{params.hyphy_dir}/runFUBAR.bf \
      $(pwd)/no_stops $(pwd)/dates $(pwd)/mrca $(pwd) $(pwd)/rates.json
    '''
}

// def make_coordinates_json(infile, outfile):
//     # FIXME: degap MRCA before running?
//     # FIXME: split up so mafft can run remotely
//     # FIXME: run mafft with --add
//     combined = os.path.join(pipeline_dir, 'combined.fasta')
//     aligned = os.path.join(pipeline_dir, 'combined.aligned.fasta')
//     cat_files([infile, globals_.config.get('Parameters', 'reference_protein')], combined)
//     mafft(combined, aligned)
//     pairwise_mrca, pairwise_ref = list(SeqIO.parse(aligned, 'fasta'))
//     ref_coords = open(globals_.config.get('Parameters', 'reference_coordinates')).read().strip().split()

//     # transfer coordinates to MRCA
//     pairwise_coords = extend_coordinates(ref_coords, str(pairwise_ref.seq))
//     mrca_coords = list(coord for char, coord in zip(str(pairwise_mrca.seq),
//                                                     pairwise_coords)
//                        if char != "-")

//     # extend MRCA coords to full MSA
//     mrca = read_single_record(infile, 'fasta', True)
//     result = extend_coordinates(mrca_coords, str(mrca.seq))

//     rdict = {'coordinates': result}
//     with open(outfile, 'w') as handle:
//         json.dump(rdict, handle, separators=(",\\n", ":"))

// def make_sequences_json(infiles, outfile):
//     alignment_file, mrca_file, coords_file = infiles
//     result = {}
//     observed = {}

//     # add sequences
//     # TODO: check if MRCA differs from our inferred MRCA
//     for r in SeqIO.parse(alignment_file, "fasta"):
//         if r.name == "Node0":
//             r.name = 'MRCA'
//         if 'ancestor' in r.name or r.name == 'MRCA':
//             if 'Ancestors' not in result:
//                 result['Ancestors'] = {}
//             result['Ancestors'][r.id] = str(r.seq)
//         else:
//             date = name_to_date(r.name)
//             if date not in observed:
//                 observed[date] = {}
//             observed[date][r.id] = str(r.seq)
//     result['Observed'] = observed

//     # add MRCA
//     mrca = read_single_record(mrca_file, 'fasta', True)
//     result['MRCA'] = str(mrca.seq)

//     # add reference
//     reference = read_single_record(globals_.config.get('Parameters', 'reference_protein'), 'fasta', True)
//     rstr = str(reference.seq)
//     with open(coords_file) as handle:
//         alignment_coords = json.load(handle)['coordinates']
//     ref_coords = open(globals_.config.get('Parameters', 'reference_coordinates')).read().strip().split()
//     cmap = dict((c, i) for i, c in enumerate(ref_coords))
//     ref_result = ''.join(list(rstr[cmap[c]] for c in alignment_coords))
//     result['Reference'] = ref_result

//     with open(outfile, 'w') as handle:
//         json.dump(result, handle, separators=(",\\n", ":"))


/* ************************************************************************** */
/* DIAGNOSIS SUB-PIPELINE */

// TODO: make this optional

/*
process diagnose {
    input:
    file "qcs*" from qcs_final_3.map{ it[0] }.collect()
    file 'hqcs' from merged_hqcs_out    

    output:
    file 'pairfile' into hqcs_qcs_pairs_out

    """
    # cat qcs
    cat qcs* > all_qcs
    
    # convert to fasta
    ${params.python} ${params.script_dir}/filter_fastx.py \
      convert fastq fasta < all_qcs > qcs_fasta

    # usearch for pairs
    ${params.usearch} -usearch_global qcs_fasta -db hqcs \
      -userout pairfile -userfields query+target \
      -id ${params.qcs_to_hqcs_identity} \
      -top_hit_only -strand both \
      -maxaccepts ${params.max_accepts} -maxrejects ${params.max_rejects} \
      -maxqt ${params.max_qcs_length_ratio} \
      -threads ${params.threads}

    # combine pairs
    mkdir alignments
    ${params.python} ${params.script_dir}/pairs_to_fasta.py \
      qcs_fasta hqcs pairfile alignments

    # align and insert gaps
    # TODO: do in parallel
    for f in alignments/*unaligned.fasta; do
        ${params.bealign} -a codon f ${f}.bam

        bam2msa ${f}.bam ${f}.bam.fasta

        ${params.python} ${params.script_dir}/insert_gaps.py \
          ${f}.bam.fasta msa ${f}.bam.fasta.gapped
    done

    # merge all
    cat alignments/gapped* > qcs_msa

    # replace gapped codons

    ${params.python} ${params.script_dir}/filter_fastx.py \
      gap_godons fasta fasta < infile > outfile
    ${params.python} ${params.script_dir}/translate.py \
      TODO

    # run diagnosis
    mkdir diagnosis
    ${params.python} ${params.script_dir}/diagnose.py hqcs qcs cn diagnosis
    """
}
*/
