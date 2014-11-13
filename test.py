import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped

from backtranslate import back_translate_gapped
from backtranslate import MissingRecord


class TestBacktranslate(unittest.TestCase):

    def compare_seqrecords(self, a, b):
        # TODO: this is probably not the best way to compare seqrecords
        self.assertEqual(str(a), str(b))
        self.assertEqual(str(a.seq.alphabet), str(b.seq.alphabet))

    def test_simple(self):
        protein = SeqRecord(Seq("MK-V", alphabet=Gapped(IUPAC.protein)))
        dna = SeqRecord(Seq("ATGAAAGTG", alphabet=IUPAC.unambiguous_dna))
        expected = SeqRecord(Seq("ATGAAA---GTG", alphabet=Gapped(IUPAC.unambiguous_dna)))
        result = back_translate_gapped(protein, dna)
        self.compare_seqrecords(expected, result)

    def test_missing_arg(self):
        protein = SeqRecord(Seq("MK-V", alphabet=Gapped(IUPAC.protein)))
        dna = SeqRecord(Seq("ATGAAAGTG", alphabet=IUPAC.unambiguous_dna))
        self.assertRaises(MissingRecord, back_translate_gapped, protein, None)
        self.assertRaises(MissingRecord, back_translate_gapped, None, dna)
        self.assertRaises(MissingRecord, back_translate_gapped, None, None)
