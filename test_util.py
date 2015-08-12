import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped

from util import insert_gaps


class TestInsertGaps(unittest.TestCase):

    def test_long_src_gap(self):
        src = "a--bc"
        target = "def"
        self.assertRaises(ValueError, insert_gaps, src, target,
                          src_gap="---", target_gap="-")

    def test_diff_lengths(self):
        src = "a-bcd"
        target = "def"
        self.assertRaises(ValueError, insert_gaps, src, target,
                          src_gap="-", target_gap="-")

