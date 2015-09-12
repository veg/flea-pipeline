#!/usr/bin/env python

"""
Takes all sequences with *perfect* reading frame, and writes the
reading frame to <infile>.perfect.fasta.

Usage:
  perfectORF [options] <infile> <outfile>
  perfectORF -h | --help

Options:
  -t INT --table=INT       Codon table [default: 1]
  -m INT --min-length=INT  Minimum protein length [default: 750]
  -v --verbose             Print progress to STDERR
  -h --help                Show this screen.

"""
import sys

from docopt import docopt

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def perfect_record(record, table, min_pro_len, verbose):
    """Returns the open reading frame, or None if no ORF found."""
    nuc = record.seq.ungap("-")
    new_str = str(nuc).upper().replace("X", "N")
    nuc = Seq(new_str, alphabet=nuc.alphabet)
    if verbose:
        sys.stderr.write(record.name)
        sys.stderr.write("\n")
        sys.stderr.write("{}...{}\n".format(nuc[:15], nuc[-15:]))
    for frame in range(3):
        # trim last few nucleotides so length is multiple of 3
        nuc_frame = nuc[frame:]
        new_len = len(nuc_frame) - (len(nuc_frame) % 3)
        nuc_frame = nuc_frame[:new_len]
        trans = nuc_frame.translate(table)
        for pro in trans.split("*"):
            if len(pro) >= min_pro_len:
                pos = trans.find(pro)
                pro_start = pro.find("M")
                if pro_start > -1:
                    if verbose:
                        sys.stderr.write(
                            "{}...{} - length {}, frame {}\n".format(
                                pro[pro_start:30], pro[-30:], len(pro), frame))
                    nuc_start = (pos * 3) + frame
                    start = nuc_start + (pro_start * 3)
                    stop = nuc_start + (len(pro) * 3)
                    orf = nuc[start:stop]
                    newrecord = SeqRecord(
                        Seq(str(orf), IUPAC.unambiguous_dna),
                        id=str(record.name), description="")
                    return newrecord
    return None


def perfect_file(infile, outfile, min_len, table=1, verbose=False):
    perfected = list(perfect_record(r, table, min_len, verbose)
                 for r in SeqIO.parse(infile, 'fasta'))
    n_bad = sum(1 for r in perfected if r is None)
    if verbose:
        sys.stderr.write("Done. Sequences without ORFs: {}\n".format(n_bad))
    perfected = list(r for r in perfected if r is not None)
    SeqIO.write(perfected, outfile, "fasta")


if __name__ == "__main__":
    args = docopt(__doc__)
    perfect_file(args["<infile>"],
                 args['<outfile'],
                 int(args["--min-length"]),
                 int(args["--table"]),
                 args["--verbose"])
