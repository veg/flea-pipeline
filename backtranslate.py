#!/usr/bin/env python
"""
Back translate protein to DNA.

Assumes sequences have the same id in each input file.

Usage:
  backtranslate.py [options] <aligned_protein> <dna> <outfile>
  backtranslate.py -h | --help

Options:
  --in-order    Match sequences by order, instead of by name.
  -h --help     Show this screen.

"""

import sys

from docopt import docopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped

from util import insert_gaps


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


def backtranslate(protein_filename, dna_filename, outfile, inorder=False):
    protein_records = SeqIO.parse(protein_filename, "fasta",
                                  alphabet=Gapped(IUPAC.protein))
    dna_records = SeqIO.parse(dna_filename, "fasta",
                              alphabet=Gapped(IUPAC.unambiguous_dna))
    if inorder:
        result_iter = (back_translate_gapped(p, d)
                       for p, d in zip(protein_records, dna_records))
    else:
        protein_dict = dict((s.id, s) for s in protein_records)
        dna_dict = dict((s.id, s) for s in dna_records)
        protein_ids = set(protein_dict)
        dna_ids = set(dna_dict)
        shared = protein_ids & dna_ids
        missing_protein = dna_ids - protein_ids
        missing_dna = protein_ids - dna_ids
        if len(missing_protein) > 0:
            sys.stderr.write('Warning: {} protein sequences'
                             ' missing.'.format(len(missing_protein)))
        if len(missing_dna) > 0:
            raise Exception('{} protein sequences have no corresponding'
                            ' dna sequence'.format(len(missing_dna)))
        result_iter = (back_translate_gapped(protein_dict[_id], dna_dict[_id])
                       for _id in shared)
    SeqIO.write(result_iter, outfile, "fasta")


if __name__ == "__main__":
    args = docopt(__doc__)
    protein_filename = args["<aligned_protein>"]
    dna_filename = args["<dna>"]
    outfile = args["<outfile>"]
    inorder = args["--in-order"]
    backtranslate(protein_filename, dna_filename, outfile, inorder=inorder)
