#!/usr/bin/env python
"""
Aligns each read to its best-matched clean read. Prints MSA to STDOUT.

Dependencies:
  - usearch
  - balign

Usage:
  align_to_clean.py <clean_gapped> <clean_ungapped> <clean_db> <dirty>
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


def usearch_global(readfile, dbfile, outdir):
    """Run usearch_global against database.

    Returns: list of (dirty id, clean id) tuples.

    """
    outfile = os.path.join(outdir, "usearch_userout.tsv")
    kwargs = dict(readfile=readfile, dbfile=dbfile, outfile=outfile)
    cmd = ("usearch -usearch_global {readfile} -db {dbfile} -id 0.8"
           " -userout {outfile} -userfields query+target -strand both")
    cmd = cmd.format(**kwargs)
    subprocess.check_call(shlex.split(cmd))
    with open(outfile) as f:
        result = list(line.strip().split("\t") for line in f.readlines())
    return result


def bealign(clean_record, dirty_records, outdir):
    """Writes records and runs bealign.

    Converts result to FASTA and return list of records.

    """
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
    return list(SeqIO.parse(outfile, "fasta"))


def align_all(dirty_filename, ungapped_filename, pairs, outdir):
    """Align dirty reads to clean reads.

    For each clean read, assembles all dirty reads matched to it and
    runs BeAlign.

    """
    match_dict = defaultdict(list)
    for dirty, clean in pairs:
        match_dict[clean].append(dirty)
    clean_records = list(SeqIO.parse(ungapped_filename, "fasta"))
    dirty_records = list(SeqIO.parse(dirty_filename, "fasta"))
    clean_dict = {r.id : r for r in clean_records}
    dirty_dict = {r.id : r for r in dirty_records}

    result = {}
    for clean_id, dirty_ids in match_dict.items():
        assert len(dirty_ids) == len(set(dirty_ids))
        target = clean_dict[clean_id]
        to_align = list(dirty_dict[i] for i in dirty_ids)
        assert len(target) % 3 == 0
        result[clean_id] = bealign(target, to_align, outdir)
    return result


def add_gaps(alignments, gapped_filename):
    pass


if __name__ == "__main__":
    args = docopt(__doc__)
    clean_gapped = args["<clean_gapped>"]
    clean_ungapped = args["<clean_ungapped>"]
    dbfile = args["<clean_db>"]
    dirty_filename = args["<dirty_fasta>"]
    outdir = "/tmp/align_{}".format(uuid.uuid4())
    os.mkdir(outdir)
    pairs = usearch_global(dirty_filename, dbfile, outdir)
    alignments = align_all(dirty_filename, clean_ungapped, pairs, outdir)
    import pdb; pdb.set_trace()
    msa = add_gaps(alignments)
    # TODO: now add gaps
    SeqIO.write(msa, sys.stdout, "fasta")
