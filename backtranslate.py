#!/usr/bin/env python
"""
Back translate protein to DNA. Prints results to STDOUT.

Assumes sequences are in order in each file.

Usage:
  backtranslate.py <gapped_protein_fasta> <ungapped_dna_fasta>
  backtranslate.py -h | --help

Options:
  -h --help     Show this screen.

"""

import sys
from itertools import zip_longest
from warnings import warn

from docopt import docopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped

from util import insert_gaps


class MissingRecord(Exception):
    pass


def back_translate_gapped(protein, dna):
    """Insert gaps from `protein` into ungapped back-translated `dna`.

    Params
    ------
    protein: SeqRecord
        Protein, with gaps, to be back-translated
    dna: SeqRecord
        DNA, without gaps, that translates to ungapped `protein`

    Returns: SeqRecord

    """
    if protein is None or dna is None:
        raise MissingRecord("One or more arguments is missing")
    if protein.id != dna.id:
        warn("Protein id={} but DNA id={}".format(protein.id, dna.id))
    try:
        dna_str = str(dna.seq.ungap())
    except ValueError:
        dna_str = str(dna.seq)
    gap_char = protein.seq.alphabet.gap_char
    gap_codon = Gapped(dna.seq.alphabet).gap_char * 3
    result_str = insert_gaps(str(protein.seq), dna_str, gap_char, gap_codon, skip=3)
    result = dna[:]
    result.seq = Seq(result_str, alphabet=Gapped(dna.seq.alphabet))
    return result


if __name__ == "__main__":
    args = docopt(__doc__)
    translated_filename = args["<gapped_protein_fasta>"]
    gapless_filename = args["<ungapped_dna_fasta>"]
    try:
        translated_handle = open(translated_filename, "rU")
        gapless_handle = open(gapless_filename, "rU")
        translated_records = SeqIO.parse(translated_handle, "fasta", alphabet=Gapped(IUPAC.protein))
        gapless_records = SeqIO.parse(gapless_handle, "fasta", alphabet=Gapped(IUPAC.unambiguous_dna))
        result_iter = (back_translate_gapped(a, b) for a, b in zip_longest(translated_records, gapless_records))
        SeqIO.write(result_iter, sys.stdout, "fasta")
    except MissingRecord:
        # manually delete generators to avoid extra exception messages
        del translated_records
        del gapless_records
        sys.stderr.write("Input files contain different numbers of sequences.\n")
        sys.exit(1)
    finally:
        translated_handle.close()
        gapless_handle.close()
