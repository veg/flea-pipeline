Dependencies
------------
- usearch
- mafft
- BioExt's bealign
- others?

Installation
------------
Still need to write a Makefile.

For now, make scripts available on the PATH.

Also, change path to references in strandCaller. (TODO: make this more
portable).

TODO: should we add references to repo, or include script to download
them?

Usage
-----

For now, we assume this is running on silverback. All the python
dependencies are already installed in a virtual environment. Just
activate it like so:

`source /opt/talign/env/bin/activate`

Run the whole pipeline with `env_pipeline.py <file>`.

The `<file>` argument is a file containing a list of fasta files, and
their sequence ids.