#!/usr/bin/env python
"""
Aligns each read to its best match in the usearch database. Prints MSA
to STDOUT.

Dependencies:
  - usearch
  - balign

Usage:
  align_to_refs.py <refs_gapped> <refs_ungapped> <refs_db> <reads>
  align_to_refs.py -h | --help

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
from Bio.Alphabet import Gapped

from util import insert_gaps


def usearch_global(readfile, dbfile, outdir):
    """Run usearch_global against database.

    Returns: list of (read id, ref id) tuples.

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


def bealign(ref_record, read_records, outdir):
    """Writes records and runs bealign.

    Converts result to FASTA and return list of records.

    """
    infile = "bealign_input_{}.fasta".format(ref_record.id)
    bamfile = "bealign_output_{}.bam".format(ref_record.id)
    outfile = "bealign_output_{}.fasta".format(ref_record.id)
    infile = os.path.join(outdir, infile)
    bamfile = os.path.join(outdir, bamfile)
    outfile = os.path.join(outdir, outfile)
    records = [ref_record]
    records.extend(read_records)
    SeqIO.write(records, infile, "fasta")

    cmd = "bealign {} {}".format(infile, bamfile)
    subprocess.check_call(shlex.split(cmd))

    cmd = "bam2msa {} {}".format(bamfile, outfile)
    subprocess.check_call(shlex.split(cmd))
    return list(SeqIO.parse(outfile, "fasta"))


def align_all(reads_filename, references_filename, pairs, outdir):
    """Align reads to ungapped clean reference reads.

    For each reference, assembles all reads matched to it and
    runs BeAlign.

    """
    match_dict = defaultdict(list)
    for read, ref in pairs:
        match_dict[ref].append(read)
    references_records = list(SeqIO.parse(references_filename, "fasta"))
    read_records = list(SeqIO.parse(reads_filename, "fasta"))
    references_dict = {r.id : r for r in references_records}
    read_dict = {r.id : r for r in read_records}

    result = []
    for ref_id, read_ids in match_dict.items():
        target = references_dict[ref_id]
        to_align = list(read_dict[i] for i in read_ids)
        result.append(bealign(target, to_align, outdir))
    return result


def add_record_gaps(reference, record):
    """Insert gaps from `reference` into another record."""
    result_str = insert_gaps(str(reference.seq), str(record.seq), "-", "-", skip=1)
    result = record[:]
    result.seq = Seq(result_str, alphabet=Gapped(record.seq.alphabet))
    return result


def add_all_gaps(alignments, gapped_filename):
    """Make local alignments match clean MSA.

    Each element of `alignments` is a list of records; the first
    record is the reference, and the rest are the reads aligned to it.

    """
    gap_dict = {r.id : r for r in SeqIO.parse(gapped_filename, 'fasta')}
    result = []
    for ref, *reads in alignments:
        gap_ref = gap_dict[ref.id]
        new_reads = list(add_record_gaps(gap_ref, r) for r in reads)
        result.extend(new_reads)
    return result


if __name__ == "__main__":
    args = docopt(__doc__)
    refs_gapped = args["<refs_gapped>"]
    refs_ungapped = args["<refs_ungapped>"]
    dbfile = args["<refs_db>"]
    reads_filename = args["<reads>"]
    outdir = "/tmp/align_{}".format(uuid.uuid4())
    os.mkdir(outdir)
    pairs = usearch_global(reads_filename, dbfile, outdir)
    alignments = align_all(reads_filename, refs_ungapped, pairs, outdir)
    msa = add_all_gaps(alignments, refs_gapped)
    lens = list(len(r) for r in msa)
    assert all(i == lens[0] for i in lens)
    SeqIO.write(msa, sys.stdout, "fasta")
