import unittest
from collections import Counter
import tempfile
import os
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped

from DNAcons import _column_consensus
from DNAcons import consensus
from DNAcons import consfile


class TestConsensus(unittest.TestCase):

    def _col_helper(self, seq, expected):
        c = Counter(seq)
        result = _column_consensus(c, seed=0)
        self.assertEquals(expected, result)

    def test_column_consensus(self):
        tests = [
            ('aabc', ('a', (0.5, 'a'))),
            ('aabbc', ('b', (2/5, 'ab'))),
            ('aabbcc', ('b', (1/3, 'abc'))),
            ('aaa', ('a', (1.0, 'a'))),
            ('a', ('a', (1.0, 'a'))),
        ]
        for seq, expected in tests:
            self._col_helper(seq, expected)

    def test_simple_consensus(self):
        result, ambi = consensus(['aab', 'abb'], seed=0)
        self.assertEquals(result, 'abb')
        self.assertEquals(ambi, ((1.0, 'a'), (0.5, 'ab'), (1.0, 'b')))

    def test_single_consensus(self):
        result, ambi = consensus(['aab'], seed=0)
        self.assertEquals(result, 'aab')
        self.assertEquals(ambi, ((1.0, 'a'), (1.0, 'a'), (1.0, 'b')))

    def test_weighted_consensus(self):
        result, ambi = consensus(['aaa', 'aaa', 'bbb'], copies=[1, 1, 8], seed=0)
        self.assertEquals(result, 'bbb')
        self.assertEquals(ambi, ((0.8, 'b'), (0.8, 'b'), (0.8, 'b')))

    def consfile_helper(self, seqs, expected, expected_ambi, ungap):
        records = list(SeqRecord(Seq(s, alphabet=Gapped(IUPAC.unambiguous_dna)))
                       for s in seqs)
        wdir = tempfile.mkdtemp()
        try:
            infile = os.path.join(wdir, 'seqs.fasta')
            outfile = os.path.join(wdir, 'cons.fasta')
            ambifile = os.path.join(wdir, 'cons.info')
            with open(infile, 'w') as handle:
                SeqIO.write(records, infile, 'fasta')
            consfile(infile, outfile, ambifile, ungap=ungap, seed=0)
            with open(outfile) as handle:
                records = list(SeqIO.parse(outfile, 'fasta'))
            self.assertEquals(len(records), 1)
            result = str(records[0].seq)
            self.assertEquals(result, expected)
            with open(ambifile) as handle:
                ambi = handle.read()
                self.assertEquals(ambi, expected_ambi)
        finally:
            shutil.rmtree(wdir)

    def test_consfile_ungapped(self):
        seqs = ["----", "ac-g", "ac-g"]
        expected = 'acg'
        exp_ambi = ''
        self.consfile_helper(seqs, expected, exp_ambi, ungap=True)

    def test_consfile_gapped(self):
        seqs = ["----", "ac-g", "ac-g"]
        expected = 'ac-g'
        exp_ambi = ''
        self.consfile_helper(seqs, expected, exp_ambi, ungap=False)

    def test_consfile_ambi(self):
        seqs = ["ac", "ag"]
        expected = 'ag'
        exp_ambi = '1 0.5 cg\n'
        self.consfile_helper(seqs, expected, exp_ambi, ungap=True)
