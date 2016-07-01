#!/usr/bin/env python
"""
Back translate protein to DNA.

Assumes sequences have the same id in each input file.

Usage:
  backtranslate.py [options] <aligned_protein> <dna> <outfile>
  backtranslate.py -h | --help

Options:
  --in-order      Match sequences by order, instead of by name.
  --remove-stops  Remove all stop codons in DNA
  -h --help       Show this screen.

"""

import sys
from flea_pipeline.util import grouper

from docopt import docopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped

from flea_pipeline.util import insert_gaps


def preprocess(protein, dna, remove_stops):
    try:
        dna = dna.ungap()
    except ValueError:
        pass
    # strip trailing nucleotides
    extra = len(dna) % 3
    if extra:
        dna = dna[:-extra]
    # remove DNA stop codon
    last_aa = dna[-3:].translate()
    if last_aa == '*':
        dna = dna[:-3]
    translated = dna.translate()
    # remove protein stop codon
    if protein[-1] == '*' and translated == str(protein.ungap())[:-1]:
        protein = protein[:-1]
    # handle internal stop codons
    if '*' in str(protein.ungap()):
        if remove_stops:
            protein = Seq(''.join(str(protein).split('*')), alphabet=protein.alphabet)
        else:
            raise Exception('protein sequence contains stop codons')
    if '*' in translated:
        if remove_stops:
            stop_positions = set(i for i, c in enumerate(translated) if c == '*')
            dna_str = ''.join(list(''.join(codon) for i, codon in enumerate(grouper(str(dna), 3))
                                   if i not in stop_positions))
            dna = Seq(dna_str, alphabet=dna.alphabet)
        else:
            raise Exception('dna sequence contains stop codons')
    if not str(dna.translate()) == str(protein.ungap()):
        raise Exception('translated sequence does not match protein')
    return protein, dna


def back_translate_gapped(protein_record, dna_record, remove_stops=False):
    """Insert gaps from `protein` into ungapped back-translated `dna`.

    Params
    ------
    protein: SeqRecord
        Protein, with gaps, to be back-translated
    dna: SeqRecord
        DNA, without gaps, that translates to ungapped `protein`

    Returns: SeqRecord

    """
    protein, dna = preprocess(protein_record.seq, dna_record.seq, remove_stops)
    gap_char = protein.alphabet.gap_char
    gap_codon = Gapped(dna.alphabet).gap_char * 3
    result_str = insert_gaps(str(protein), str(dna), gap_char, gap_codon)
    result = dna_record[:]
    result.seq = Seq(result_str, alphabet=Gapped(dna.alphabet))
    return result


def backtranslate(protein_filename, dna_filename, outfile,
                  inorder=False, remove_stops=False):
    protein_records = SeqIO.parse(protein_filename, "fasta",
                                  alphabet=Gapped(IUPAC.protein))
    dna_records = SeqIO.parse(dna_filename, "fasta",
                              alphabet=Gapped(IUPAC.unambiguous_dna))
    if inorder:
        result_iter = (back_translate_gapped(p, d, remove_stops)
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
        result_iter = (back_translate_gapped(protein_dict[_id],
                                             dna_dict[_id],
                                             remove_stops)
                       for _id in shared)
    SeqIO.write(result_iter, outfile, "fasta")


if __name__ == "__main__":
    args = docopt(__doc__)
    protein_filename = args["<aligned_protein>"]
    dna_filename = args["<dna>"]
    outfile = args["<outfile>"]
    inorder = args["--in-order"]
    remove_stops = args["--remove-stops"]
    backtranslate(protein_filename, dna_filename, outfile,
                  inorder=inorder, remove_stops=remove_stops)
