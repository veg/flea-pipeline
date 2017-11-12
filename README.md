Overview
========

FLEA is a bioinformatics pipeline for analyzing longitudinal
sequencing data from the Pacific Biosciences RS-II or Sequel. It
currently supports full-length HIV *env* sequences.

The pipeline takes a set of FASTQ files, one per time point,
containing circular consensus sequence (CCS) reads, which can be
obtained using the ”Reads of Insert“ protocol on PacBio’s SMRTportal
or SMRTanalysis tools. It produces a JSON file containing the
following results:

- a multiple sequence alignment of high-quality consensus sequences
  for each time point

- a maximum-likelihood phylogenetic tree, inferred using
  [FastTree](http://www.microbesonline.org/fasttree/)

- the most recent common ancestor (MRCA) and other inferred ancestor
  sequences

- a two-dimensional embedding that respects TN93 sequence distances

- per-site selection pressure, inferred using
  [FUBAR](https://veg.github.io/hyphy-site/methods/selection-methods/),
  and other per-site evolutionary metrics

- per-segment evolutionary and phenotypic metrics, inferred using
  [HyPhy](http://www.hyphy.org/)

The pipeline logic is implemented in
[Nextflow](https://www.nextflow.io/). A full description of the
pipeline has been submitted for publication. A link to the journal
article will be added here when it is available.

Setup
=====

Dependencies
------------
- [Nextflow](https://www.nextflow.io/)
- [Python](https://www.python.org/)
- [USEARCH](https://www.drive5.com/usearch/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [HyPhy](http://www.hyphy.org/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [TN93](https://github.com/veg/tn93)
- [GNU parallel](https://www.gnu.org/software/parallel/)
- Python dependencies (see below)

Install Python scripts
----------------------

FLEA comes with a `flea` Python package containing scripts used
throughout the pipeline. To install requirements and the scripts
themselves (virtualenv recommended):

    pip install -r requirements.txt
    pip install -r requirements2.txt
    python setup.py install

To test:

    python setup.py nosetests


Configuration
-------------

The default config file is `nextflow.config`. It is recommended that
you make a seperate config file that overrides any options that need
to be customized. For more information on Nextflow-specific
configuration, see [the Nextflow
documentation](https://www.nextflow.io/docs/latest/config.html).

At the very least, `params.reference_dir` and the parameters that
depend on it need to point to the various reference files used by the
pipeline:

- `params.reference_db`: FASTA file of reference sequences
- `params.contaminants_db`: FASTA file of contaminant sequences
- `params.reference_dna`: reference DNA sequence
- `params.reference_protein`: reference amino acid sequence
- `params.reference_coordinates`: 


Usage
=====

Write a control file containing a list of FASTQ files, visit codes,
and dates, seperated by spaces.

    <file> <visit code> <date>
    <file> <visit code> <date>
    ....

Dates must be in 'YYYYMMDD' format.

Run the pipeline with Nextflow:

    nextflow path/to/flea.nf -c path/to/custom/config/file \
      --infile path/to/metadata \
      --results_dir path/to/results

The results directory will contain output from lots of pipeline
steps. The two files that contain the final results are:

- `session.json`: a JSON file to be visualized with
  [`flea-web-app`](https://github.com/veg/flea-web-app).

- `session.zip`: a zip file with FASTA files for the consensus
  sequences, ancestors, and MRCA, and a Newick file containing the
  rooted phylogenetic tree.