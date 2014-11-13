#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped


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
    gap_char = protein.seq.alphabet.gap_char
    gap_codon = Gapped(dna.seq.alphabet).gap_char * 3
    try:
        dna_str = str(dna.seq.ungap())
    except ValueError:
        dna_str = str(dna.seq)
    codons = []
    for a in protein.seq:
        if a == gap_char:
            codons.append(gap_codon)
        else:
            codons.append(dna_str[:3])
            dna_str = dna_str[3:]
    result = dna[:]
    result.seq = Seq("".join(codons), alphabet=Gapped(dna.seq.alphabet))
    return result


if __name__ == "__main__":
    import sys

    args = sys.argv[1:]
    if len(args) != 2:
        sys.stderr.write("usage: backtranslate.py <gapped protein fasta> <ungapped dna fasta>\n")
        sys.stderr.write("       prints result to STDOUT\n")
        sys.exit()

    translated_filename = args[0]
    gapless_filename = args[1]

    translated_handle = open(translated_filename, "rU")
    gapless_handle = open(gapless_filename, "rU")
    translated_records = SeqIO.parse(translated_handle, "fasta", alphabet=Gapped(IUPAC.protein))
    gapless_records = SeqIO.parse(gapless_handle, "fasta", alphabet=Gapped(IUPAC.unambiguous_dna))
    result_iter = (back_translate_gapped(a, b) for a, b in zip(translated_records, gapless_records))
    SeqIO.write(result_iter, sys.stdout, "fasta")
    translated_handle.close()
    gapless_handle.close()
