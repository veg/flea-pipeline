import unittest

from correct_shifts import correct_shifts

class TestCorrectShifts(unittest.TestCase):

    def _test(self, seq, ref, expected, keep=False):
        self.assertEquals(correct_shifts(seq, ref, keep=keep), expected)

    def test_correct_shifts(self):
        cases = (
            # deletions
            ('AC---T', 'ACGACG', 'ACT'),

            # insertions
            ('ACGT', 'AC-T', 'ACT'),
            ('ACGTT', 'AC--T', 'ACT'),
            ('ACGTGT', 'AC---T', 'ACGTGT'),
            ('ACAAAG', 'AC---G', 'ACAAAG'),
            ('ACAAAAG', 'AC----G', ''),
            ('ACAAACCCG', 'AC------G', 'ACAAACCCG'),

            # both; keep=False
            ('AG-GTTT', 'ACC-TTT', ''),
        )
        for seq, ref, expected in cases:
            self._test(seq, ref, expected)

    def test_keep(self):
        cases = (
            ('ACAAAAG', 'AC----G', 'ACAAAAG'),

            # deletions
            ('AC-ACG', 'ACGACG', 'ACXACG'),
            ('AC--CG', 'ACGACG', 'ACXXCG'),
            ('AC---T', 'ACGACG', 'ACT'),
            ('AC----ACG', 'ACGACGACG', 'ACXACG'),

            # both; keep=True
            ('AC-TTTT', 'ACCTTT-', 'ACXTTT'),
        )
        for seq, ref, expected in cases:
            self._test(seq, ref, expected, keep=True)


    def test_mismatched_lengths(self):
        seq = 'AC'
        ref = 'ACC'
        self.assertRaises(ValueError, correct_shifts, seq, ref)

    def test_reference_length(self):
        seq = 'AC'
        ref = 'ACCC'
        self.assertRaises(ValueError, correct_shifts, seq, ref)


if __name__ == '__main__':
    unittest.main()
