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

params.infile = "$HOME/flea/data/P018/data/P018_3mROIx6.fastq"
params.time_label = "V03"

inputfile = Channel.fromPath(params.infile)

process quality_filter {
    input:
    file 'ccs.fastq' from inputfile

    output:
    file qfiltered

    """
    ${params.usearch} -fastq_filter ccs.fastq \
      -fastq_maxee_rate ${params.max_error_rate} \
      -threads ${params.threads} \
      -fastq_qmax ${params.qmax} -fastq_minlen ${params.min_len} \
      -fastqout qfiltered \
      -relabel "${params.time_label}_ccs_"
    """
}

process trim_heads_and_tails {
    input:
    file qfiltered

    output:
    file trimmed

    "${params.python} ${params.script_dir}/trim.py --jobs ${params.threads} --fastq qfiltered trimmed"
}

process filter_runs {
    input:
    file trimmed

    output:
    file no_runs

    """
    #!/usr/bin/env python

    from Bio import SeqIO
    from flea.util import run_regexp

    cregexp = run_regexp(${params.run_length})
    records = SeqIO.parse("trimmed", 'fastq')
    result = (r for r in records if cregexp.search(str(r.seq)) is None)
    SeqIO.write(result, "no_runs", 'fastq')
    """
}

process filter_contaminants {
    input:
    file no_runs

    output:
    file uncontam

    """
    ${params.usearch} -usearch_global no_runs \
      -db ${params.contaminants_db} -id ${params.contaminant_identity} -strand both \
      -threads ${params.threads} -notmatchedfq uncontam
    """
}

process filter_refdb {
    input:
    file uncontam

    output:
    file matches_db
    file userout

    """
    ${params.usearch} -usearch_global uncontam -db ${params.reference_db} \
      -id ${params.reference_identity} \
      -userfields qstrand+tstrand+caln \
      -top_hit_only -strand both \
      -maxaccepts ${params.max_accepts} -maxrejects ${params.max_rejects} \
      -threads ${params.threads} \
      -matchedfq matches_db -userout userout
    """
}

process trim_terminal_gaps {
    input:
    file matches_db
    file userout

    output:
    file gaptrimmed

     "${params.python} ${params.script_dir}/trim_terminal_gaps.py matches_db userout gaptrimmed"
}

// FIXME: do not hardcode min/max length
process filter_length {
    input:
    file gaptrimmed

    output:
    file lenfiltered

    """
    #!/usr/bin/env python

    from Bio import SeqIO

    records = SeqIO.parse("gaptrimmed", 'fastq')
    result = (r for r in records if ${params.min_len} <= len(r.seq) <= ${params.max_len})
    SeqIO.write(result, "lenfiltered", 'fastq')

    """
}
