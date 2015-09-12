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

from flea_pipeline.DNAcons import _column_consensus
from flea_pipeline.DNAcons import consensus
from flea_pipeline.DNAcons import consfile


class TestConsensus(unittest.TestCase):

    def _col_helper(self, seq, expected):
        c = Counter(seq)
        result = _column_consensus(c, seed=0)
        self.assertEquals(expected, result)

    def test_column_consensus(self):
        tests = [
            ('aabc', ('a', (0.5, ['a']))),
            ('aabbc', ('b', (2/5, ['a', 'b']))),
            ('aabbcc', ('b', (1/3, ['a', 'b', 'c']))),
            ('aaa', ('a', (1.0, ['a']))),
            ('a', ('a', (1.0, ['a']))),
        ]
        for seq, expected in tests:
            self._col_helper(seq, expected)

    def test_simple_consensus(self):
        result, ambi = consensus(['aab', 'abb'], seed=0)
        self.assertEquals(result, 'abb')
        self.assertEquals(ambi, ((1.0, ['a']), (0.5, ['a', 'b']), (1.0, ['b'])))

    def test_single_consensus(self):
        result, ambi = consensus(['aab'], seed=0)
        self.assertEquals(result, 'aab')
        self.assertEquals(ambi, ((1.0, ['a']), (1.0, ['a']), (1.0, ['b'])))

    def test_weighted_consensus(self):
        result, ambi = consensus(['aaa', 'aaa', 'bbb'], copies=[1, 1, 8], seed=0)
        self.assertEquals(result, 'bbb')
        self.assertEquals(ambi, ((0.8, ['b']), (0.8, ['b']), (0.8, ['b'])))

    def test_weighted_non_codon(self):
        result, ambi = consensus(['aac', 'acc', 'bbb'], copies=[3, 3, 4], seed=0)
        self.assertEquals(result, 'abc')
        self.assertEquals(ambi, ((0.6, ['a']), (0.4, ['b']), (0.6, ['c'])))

    def test_weighted_codon(self):
        result, ambi = consensus(['aacddd', 'accddd', 'bbbddd'], copies=[3, 3, 4], codon=True, seed=0)
        self.assertEquals(result, 'bbbddd')
        self.assertEquals(ambi, ((0.4, ['bbb']), (1.0, ['ddd'])))

    def consfile_helper(self, seqs, expected, copynumbers, expected_ambi, ungap, codon=False):
        records = list(SeqRecord(Seq(s, alphabet=Gapped(IUPAC.unambiguous_dna)),
                                 id="seq_{}".format(i))
                       for i, s in enumerate(seqs))
        cdict = dict((r.id, n) for r, n in zip(records, copynumbers))
        wdir = tempfile.mkdtemp()
        try:
            infile = os.path.join(wdir, 'seqs.fasta')
            copynumber_file = os.path.join(wdir, 'copynumbers')
            outfile = os.path.join(wdir, 'cons.fasta')
            ambifile = os.path.join(wdir, 'cons.info')
            SeqIO.write(records, infile, 'fasta')
            SeqIO.write(records, infile, 'fasta')
            with open(copynumber_file, 'w') as handle:
                for k, val in cdict.items():
                    handle.write("{}\t{}\n".format(k, val))

            consfile(infile, outfile, ambifile=ambifile,
                     copynumber_file=copynumber_file,
                     ungap=ungap, codon=codon, seed=0)

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
        copynumbers = [1, 1, 1]
        exp_ambi = ''
        self.consfile_helper(seqs, expected, copynumbers, exp_ambi, ungap=True)

    def test_consfile_gapped(self):
        seqs = ["----", "ac-g", "ac-g"]
        expected = 'ac-g'
        copynumbers = [1, 1, 1]
        exp_ambi = ''
        self.consfile_helper(seqs, expected, copynumbers, exp_ambi, ungap=False)

    def test_consfile_ambi(self):
        seqs = ["ac", "ag"]
        expected = 'ag'
        copynumbers = [1, 1]
        exp_ambi = "1 0.5 ['c', 'g']\n"
        self.consfile_helper(seqs, expected, copynumbers, exp_ambi, ungap=True)

    def test_consfile_codon(self):
        seqs = ["aac", "acc", "bbb"]
        expected = 'bbb'
        copynumbers = [2, 2, 3]
        exp_ambi = ""
        self.consfile_helper(seqs, expected, copynumbers, exp_ambi, ungap=True, codon=True)
