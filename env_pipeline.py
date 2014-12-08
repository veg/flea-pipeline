#!/usr/bin/env python

"""
Runs the complete pipeline on a set of fasta files.

Input is a file with the following format:

<file_1> <id_1>
<file_2> <id_2>
...
<file_n> <id_n>


Usage:
  env_pipeline [options] <file>
  env_pipeline -h | --help

Options:
  --percent-identity=FLOAT  Percent identity for clustering [default: 0.99]
  --discard-lb=INT          Lower bound for selecting clusters [default: 5]
  --subsample-ub=INT        Upper bound for selecting clusters [default: 30]
  -v --verbose              Print progress.
  -h --help                 Show this screen.

"""

from subprocess import check_call
import shlex
from os import path
import csv
from collections import namedtuple

#FIX THIS PLS
import os

from docopt import docopt

from translate import translate
from backtranslate import backtranslate
from align_to_refs import align_to_refs


def add_suffix(filename, suffix):
    """Add suffix before extension

    >>> add_suffix("path/to/filename.ext", "processed")
    "path/to/filename_processed.ext"

    """
    base, ext = path.splitext(filename)
    return "".join([base, "_", suffix, ext])


def call(cmd, **kwargs):
    check_call(shlex.split(cmd), **kwargs)


def flatten(it):
    return (b for a in it for b in a)


def cat_files(files, outfile, chunk=2**14):
    """Concatenate multiple files in chunks."""
    out_handle = open(outfile, "w")
    for f in files:
        handle = open(f)
        while True:
            data = handle.read(chunk)
            if data:
                out_handle.write(data)
            else:
                break


Timepoint = namedtuple('Timepoint', ['file', 'id', 'date'])


if __name__ == "__main__":
    args = docopt(__doc__)

    percent_identity = args["--percent-identity"]
    discard_lb = args["--discard-lb"]
    subsample_ub = args["--subsample-ub"]
    is_verbose = args["--verbose"]

    # FIXME: for now, filenames cannot have spaces. Make this more
    # robust. For instance use tab-seperated values.
    with open(args["<file>"], newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')
        timepoints = list(Timepoint(f, i, d) for f, i, d in reader)

    # TODO: do this in parallel
    for t in timepoints:
        cmd = ("processTimestep {f} {seq_id} {percent_identity}"
               " {discard_lb} {subsample_ub}")
        call(cmd.format(f=t.file, seq_id=t.id, percent_identity=percent_identity,
                        discard_lb=discard_lb, subsample_ub=subsample_ub))

    # append all perfect orfs and make database
    perfect_files = list("{}.collapsed.fasta.perfect.fasta".format(t.id)
                         for t in timepoints)
    data_dir = path.dirname(path.abspath(timepoints[0].file))  # assume all are in same directory
    all_orfs_file = path.join(data_dir, "all_perfect_orfs.fasta")
    db_file = path.join(data_dir, "all_perfect_orfs.udb")
    cat_files(perfect_files, all_orfs_file)
    call("usearch -makeudb_usearch {} -output {}".format(all_orfs_file, db_file))

    # codon alignment of all perfect orfs
    translated_file = add_suffix(all_orfs_file, "translated")
    aligned_file = add_suffix(translated_file, "aligned")
    backtranslated_file = add_suffix(aligned_file, "backtranslated")
    translate(all_orfs_file, translated_file)
    call("mafft --auto {}".format(translated_file), stdout=open(aligned_file, "w"))
    backtranslate(aligned_file, all_orfs_file, backtranslated_file)

    # run alignment in each timestep
    for t  in timepoints:
        # TODO: fix these names
        infile = "".join([t.file, ".pbformatfixed.fastq.good.fasta.matches.fasta.seconds.fasta"])
        outfile = "".join([t.file, "_ALIGNED.fasta"])
        align_to_refs(infile, all_orfs_file, backtranslated_file,
                      db_file, outfile)
    #FIX THIS HORRIBLE NONSENSE PLS                  
    os.system("cat *ALIGNED.fasta >allTimepoints.ALIGNED.fasta")
    os.system("trees ALIGNED")
