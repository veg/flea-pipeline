Python dependencies
-------------------
- BioPython
- docopt
- ruffus


3rd party dependencies
----------------------
- usearch (version >= 8.0.1570)
- mafft
- BioExt's bealign


Setup
-----
1. Copy `env_pipeline.config.tpl` to `env_pipeline.config`.

2. Download contaminants and references files.

3. Generate usearch databases for each. (http://www.drive5.com/usearch/manual/cmd_makeudb_usearch.html)

4. Edit paths in `env_pipeline.config` to point to those databases.


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

If no config file is specified, we look for it in in the data
directory, and then in the script directory.