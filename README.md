Python dependencies
-------------------
- BioPython
- docopt
- ruffus


3rd party dependencies
----------------------
- usearch
- mafft
- BioExt's bealign


Setup
-----
1. Download contaminants and references files.

2. Generate usearch databases for each. (http://www.drive5.com/usearch/manual/cmd_makeudb_usearch.html)

3. Edit paths in `env_pipeline.config` to point to those databases.


Usage
-----
For now, we assume this is running on silverback. All the Python
dependencies are already installed in a virtual environment. Just
activate it like so:

`source /opt/talign/env/bin/activate`

Run the whole pipeline with `env_pipeline.py <control file>`.

The `<control file>` argument is a file containing a list of fasta
files, their sequence ids, and their dates, seperated by spaces.

    <file> <id> <date>
    <file> <id> <date>
    ....

Dates must be in 'YYYYMMDD' format.
