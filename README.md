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
TODO: script to download/build references.

1. Download contaminants references files.

2. Generate usearch databases for each. (http://www.drive5.com/usearch/manual/cmd_makeudb_usearch.html)

3. Edit paths in `config` to point to those databases.


Usage
-----

For now, we assume this is running on silverback. All the python
dependencies are already installed in a virtual environment. Just
activate it like so:

`source /opt/talign/env/bin/activate`

Run the whole pipeline with `env_pipeline.py <file>`.

The `<file>` argument is a file containing a list of fasta files, and
their sequence ids.