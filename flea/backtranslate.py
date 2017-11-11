#!/usr/bin/env python
"""
Back translate protein to DNA.

Assumes sequences have the same id in each input file.

Tolerate stop codons in DNA as long as matching '*' appears in
protein.

"""

import sys

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped

from flea.util import grouper
from flea.util import insert_gaps


def preprocess(protein, dna):
    try:
        dna = dna.ungap()
    except ValueError:
        pass
    # strip trailing nucleotides
    extra = len(dna) % 3
    if extra:
        dna = dna[:-extra]
    translated = dna.translate()
    # remove trailing protein stop codon
    if protein[-1] == '*' and translated == str(protein.ungap())[:-1]:
        protein = protein[:-1]
    if str(dna.translate()) != str(protein.ungap()):
        raise Exception('translated sequence does not match protein')
    return protein, dna


def back_translate_gapped(protein_record, dna_record):
    """Insert gaps from `protein` into ungapped back-translated `dna`.

    Params
    ------
    protein: SeqRecord
        Protein, with gaps, to be back-translated
    dna: SeqRecord
        DNA, without gaps, that translates to ungapped `protein`

    Returns: SeqRecord

    """
    protein, dna = preprocess(protein_record.seq, dna_record.seq)
    gap_char = protein.alphabet.gap_char
    gap_codon = Gapped(dna.alphabet).gap_char * 3
    result_str = insert_gaps(str(protein), str(dna), gap_char, gap_codon)
    result = dna_record[:]
    result.seq = Seq(result_str, alphabet=Gapped(dna.alphabet))
    return result


@click.command()
@click.argument('protein_filename')
@click.argument('dna_filename')
@click.argument('outfile')
@click.option('--in-order', is_flag=True,
              help='Match sequences by order, instead of by name.')
def backtranslate(protein_filename, dna_filename, outfile, in_order=False):
    protein_records = SeqIO.parse(protein_filename, "fasta",
                                  alphabet=Gapped(IUPAC.protein))
    dna_records = SeqIO.parse(dna_filename, "fasta",
                              alphabet=Gapped(IUPAC.unambiguous_dna))
    if in_order:
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
        result_iter = (back_translate_gapped(protein_dict[_id],
                                             dna_dict[_id])
                       for _id in sorted(shared))
    SeqIO.write(result_iter, outfile, "fasta")


if __name__ == "__main__":
    backtranslate()
