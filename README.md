Dependencies
------------
- Nextflow
- Python
- usearch
- mafft
- HyPhy
- TN93
- GNU parallel

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


Usage
-----
Write a control file containing a list of fastq files, their sequence ids, and
their dates, seperated by spaces.

    <file> <label> <date>
    <file> <label> <date>
    ....

Dates must be in 'YYYYMMDD' format.

Run the pipeline with Nextflow:

    nextflow path/to/flea.nf --infile path/to/metadata --results_dir path/to/results
