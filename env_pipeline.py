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
from subprocess import Popen, PIPE, STDOUT
import shlex
from os import path
import csv
import json
from collections import namedtuple
from datetime import datetime
import sys

#FIX THIS PLS
import os

from docopt import docopt

from Bio import SeqIO

from translate import translate
from backtranslate import backtranslate
from align_to_refs import align_to_refs
from DNAcons import dnacons


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


def mrca(infile, recordfile, outfile, oldest_id):
    """Writes records from `infile` to `recordfile` that have an ID
    corresponding to `oldest_id`. Then runs DNAcons, writing result to
    `outfile`.

    """
    len_id = len(oldest_id)
    records = SeqIO.parse(infile, "fasta")
    oldest_records = (r for r in records if r.id.startswith(oldest_id))
    SeqIO.write(oldest_records, recordfile, "fasta")
    dnacons(recordfile, id_str="mrca", outfile=outfile)


Timepoint = namedtuple('Timepoint', ['file', 'id', 'date'])


if __name__ == "__main__":
    args = docopt(__doc__)

    percent_identity = args["--percent-identity"]
    discard_lb = args["--discard-lb"]
    subsample_ub = args["--subsample-ub"]
    is_verbose = args["--verbose"]

    # FIXME: parts of this script assumes CWD is the data
    # directory. Make it not rely on that.

    # FIXME: for now, filenames cannot have spaces. Make this more
    # robust. For instance use tab-seperated values.
    with open(args["<file>"], newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')
        timepoints = list(Timepoint(f, i, d) for f, i, d in reader)

    # assume all are in same directory
    data_dir = path.dirname(path.abspath(timepoints[0].file))

    # TODO: do this in parallel
    for t in timepoints:
        cmd = ("processTimestep {f} {seq_id} {percent_identity}"
               " {discard_lb} {subsample_ub}")
        call(cmd.format(f=t.file, seq_id=t.id, percent_identity=percent_identity,
                        discard_lb=discard_lb, subsample_ub=subsample_ub))

    # append all perfect orfs and make database
    perfect_files = list("{}.collapsed.fasta.perfect.fasta".format(t.id)
                         for t in timepoints)
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

    # write the dates file for frontend - all perfect sequences
    # NOTE: we assume the first part of the record id is the timestamp
    # id, followed by a dot.
    dates_file = os.path.join(data_dir, "merged.dates")
    records = list(SeqIO.parse(backtranslated_file, "fasta"))
    id_to_date = {t.id : t.date for t in timepoints}
    with open(dates_file, "w") as handle:
        outdict = {r.id: id_to_date[r.id.split(".")[0]] for r in records}
        json.dump(outdict, handle, separators=(",\n", ":"))

    # get a DNA consensus for perfect ORF corresponding to earliest date
    strptime = lambda t: datetime.strptime(t.date, "%Y%m%d")
    oldest_timepoint = min(timepoints, key=strptime)
    oldest_records_filename = add_suffix(backtranslated_file,
                                         "oldest_{}".format(oldest_timepoint.id))
    mrca_filename = os.path.join(data_dir, "mrca.seq")
    mrca(backtranslated_file, oldest_records_filename, mrca_filename,
         oldest_timepoint.id)

    # run alignment in each timestep
    timestep_aligned_files = []
    for t  in timepoints:
        # TODO: fix these names
        infile = "".join([t.file, ".pbformatfixed.fastq.good.fasta.matches.fasta.seconds.fasta"])
        outfile = "".join([t.file, "_ALIGNED.fasta"])
        align_to_refs(infile, all_orfs_file, backtranslated_file,
                      db_file, outfile)
        timestep_aligned_files.append(outfile)

    # concatenate all aligned timesteps
    final_alignment_file = os.path.join(data_dir, "allTimepoints.ALIGNED.fasta")
    cat_files(timestep_aligned_files, final_alignment_file)

    # get HXB2 coordinates -> merged.prot.parts
    # TODO: result is saved in script directory, not data directory
    script_dir = os.path.split(__file__)[0]
    script_name = 'HXB2partsSplitter.bf'
    hxb2_script = os.path.join(script_dir, script_name)
    cmd = 'HYPHYMP {}'.format(hxb2_script)
    p = Popen(shlex.split(cmd), stdin=PIPE)
    script_input = " {}\n".format(final_alignment_file).encode()
    p.communicate(input=script_input)
    returncode = p.wait()
    if not returncode:
        sys.stderr.write("{} failed.".format(script_name))

    #FIX THIS HORRIBLE NONSENSE PLS
    os.system("trees ALIGNED")
