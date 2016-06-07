Python dependencies
-------------------
pip-tools is used to list requirements and pin versions.

Install pip-tools. Then just run 'pip-sync requirements.txt'

Then run `pip install -r requirements2.txt`. This is done seperately
because BioExt needs NumPy to be present.


3rd party dependencies
----------------------
- usearch (version >= 8.1.1861)
- mafft
- Julia
- pbs-drmaa

To install Julia dependencies:

    Pkg.clone("git@bitbucket.org:keren/Quiver2.jl.git")
    Pkg.resolve()

To install pbs-drmaa: do the usual

    configure --prefix=$HOME; make; make install

then ensure that the environment variable DRMAA_LIBRARY_PATH points to libdrmaa.so.1:

     export DRMAA_LIBRARY_PATH=$HOME/lib/libdrmaa.so.1


Installation
------------
`python setup.py install`


Usage
-----
1. Copy `flea.config.tpl` to `flea.config` in your data directory.

2. Edit other parameters in `flea.config`.

3. Run the whole pipeline with `flea.py <control file>`.

The `<control file>` argument is a file containing a list of fasta
files, their sequence ids, and their dates, seperated by spaces.

    <file> <label> <date>
    <file> <label> <date>
    ....

Dates must be in 'YYYYMMDD' format.

If no config file is specified, we look for it in in the data
directory.

Pre-existing alignment
----------------------

Control file should contain:

    <prefix> <label> <date>
    <prefix> <label> <date>
    ...

where each sequence id begins with a timepoint-unique `<prefix>`.

Run with the option `--alignment <alignment file>`.

Also, optionally add the option `--copynumbers <copynumber file>`,
where each line of the copynumbers file contains a tab-seperated
sequence id and integer:

    <sequence id>\t<integer>
    ...

Testing
-------

`python setup.py nosetests`