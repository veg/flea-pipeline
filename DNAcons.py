#!/usr/bin/env python

"""
Prints a gapless consensus sequence with any nucleotides lower that
50% identity recieving an ambiguous character X.

Sequences must be aligned, and all the same length'

Usage:
  DNAcons.py [options] <infile>
  DNAcons.py -h | --help

Options:
  -v --verbose           Print progress to STDERR
  --keep-gaps            Do not ungap the consensus
  --id=<STRING>          Record id for the fasta output
  -o --outfile=<STRING>  Name of output file
  -h --help              Show this screen

"""

import sys

from docopt import docopt

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC


def dnacons(filename, outfile=None, id_str=None, ambiguous=None, ungap=True, verbose=False):
    if ambiguous is None:
        ambiguous = 'N'
    alignment = AlignIO.read(filename, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.gap_consensus(threshold=0.00, ambiguous=ambiguous)
    if ungap:
        consensus = consensus.ungap("-")

    # use the first entry for the name, just so know what is what later
    rec = alignment[0]
    if id_str is None:
        id_str = "{}_cons".format(rec.name)
    if outfile is None:
        outfile = "{}.cons.fasta".format(filename)
    writer = FastaWriter(open(outfile, "w"), wrap=None)
    writer.write_header()
    newrecord = SeqRecord(Seq(str(consensus), IUPAC.unambiguous_dna), id=id_str, description="")
    writer.write_record(newrecord)
    if verbose:
        sys.stderr.write('Consensus written with name {}.\n'.format(id_str))


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["<infile>"]
    verbose = args["--verbose"]
    id_str = args["--id"]
    outfile = args["--outfile"]
    keep_gap = args["--keep-gaps"]
    ungap = not keep_gap
    dnacons(filename, outfile, id_str, ungap=ungap, verbose=verbose)

