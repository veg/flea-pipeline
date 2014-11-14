#!/usr/bin/env python
"""
Aligns each read to its best-matched clean read. Prints MSA to STDOUT.

Dependencies:
  - usearch
  - balign

Usage:
  align_to_clean.py <clean_fasta> <dirty_fasta>
  align_to_clean.py -h | --help

Options:
  -h --help      Show this screen.

"""
import itertools
import subprocess
import shlex
import uuid
import os
from collections import defaultdict
import sys

from docopt import docopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


# Taken from itertools recipes
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks."
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def make_db(infile, dbfile):
    """Make usearch database of reads in `infile`."""
    kwargs = dict(infile=infile, dbfile=dbfile)
    cmd = "usearch -makeudb_usearch {infile} -output {dbfile}".format(**kwargs)
    subprocess.check_call(shlex.split(cmd))


def usearch_global(readfile, dbfile, pairfile):
    """Run usearch, outputting sequence ids for matches."""
    kwargs = dict(readfile=readfile, dbfile=dbfile, pairfile=pairfile)
    cmd = ("usearch -usearch_global {readfile} -db {dbfile} -id 0.8"
           " -userout {pairfile} -userfields query+target -strand both")
    cmd = cmd.format(**kwargs)
    subprocess.check_call(shlex.split(cmd))


def match_reads(clean_filename, dirty_filename, outdir):
    """For each dirty read, find the best clean read.

    Uses usearch to build a database of clean reads and run global
    alignment.

    Returns: list of (dirty id, clean id) tuples.

    """
    dbfile = os.path.join(outdir, "clean_database.udb")
    pairfile = os.path.join(outdir, "usearch_userout.tsv")
    make_db(clean_filename, dbfile)
    usearch_global(dirty_filename, dbfile, pairfile)
    with open(pairfile) as f:
        result = list(line.strip().split("\t") for line in f.readlines())
    return result


def read_fasta(filename):
    with open(filename, "rU") as handle:
        return list(SeqIO.parse(handle, "fasta"))


def write_fasta(records, filename):
    with open(filename, "w") as handle:
        SeqIO.write(records, handle, "fasta")

def make_record_dict(*records):
    result = dict()
    for a in records:
        for r in a:
            result[r.id] = r
    return result


def align_to_clean(clean_record, dirty_records, outdir):
    """Write records and run BeAlign. Convert result to FASTA and return list of records."""
    infile = "bealign_input_{}.fasta".format(clean_record.id)
    bamfile = "bealign_output_{}.bam".format(clean_record.id)
    outfile = "bealign_output_{}.fasta".format(clean_record.id)
    infile = os.path.join(outdir, infile)
    bamfile = os.path.join(outdir, bamfile)
    outfile = os.path.join(outdir, outfile)
    records = [clean_record]
    records.extend(dirty_records)
    write_fasta(records, infile)

    cmd = "bealign {} {}".format(infile, bamfile)
    subprocess.check_call(shlex.split(cmd))

    cmd = "bam2msa {} {}".format(bamfile, outfile)
    subprocess.check_call(shlex.split(cmd))
    return read_fasta(outfile)


def align_all(clean_filename, dirty_filename, pairs, outdir):
    """Align dirty reads to clean reads.

    For each clean read, assembles all dirty reads matched to it and
    runs BeAlign.

    """
    clean_id_to_dirty_ids = defaultdict(list)
    for dirty, clean in pairs:
        clean_id_to_dirty_ids[clean].append(dirty)
    clean_records = read_fasta(clean_filename)
    dirty_records = read_fasta(dirty_filename)
    clean_record_dict = make_record_dict(clean_records)
    dirty_record_dict = make_record_dict(dirty_records)

    result = []
    for clean_id, dirty_ids in clean_id_to_dirty_ids.items():
        target = clean_record_dict[clean_id]
        to_align = list(dirty_record_dict[i] for i in dirty_ids)
        result.append(align_to_clean(target, to_align, outdir))
    return result


if __name__ == "__main__":
    args = docopt(__doc__)
    clean_filename = args["<clean_fasta>"]
    dirty_filename = args["<dirty_fasta>"]
    outdir = "/tmp/align_{}".format(uuid.uuid4())
    os.mkdir(outdir)
    pairs = match_reads(clean_filename, dirty_filename, outdir)
    msa = align_all(clean_filename, dirty_filename, pairs, outdir)
    SeqIO.write(msa, sys.stdout, "fasta")
