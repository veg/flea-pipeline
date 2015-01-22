import unittest

from correct_shifts import correct_shifts

class TestCorrectShifts(unittest.TestCase):

    def _test(self, seq, ref, expected):
        self.assertEquals(correct_shifts(seq, ref), expected)

    def test_correct_shifts(self):
        cases = (
            # deletions
            ('AC-T', 'ACGT', 'ACXT'),
            ('AC--T', 'ACGGT', 'ACXXT'),
            ('AC---T', 'ACGGGT', 'ACT'),
            ('AC----T', 'ACGGGGT', 'ACXT'),

            # insertions
            ('ACGT', 'AC-T', 'ACT'),
            ('ACGTT', 'AC--T', 'ACT'),
            ('ACGTGT', 'AC---T', 'ACGTGT'),
            ('ACAAAGT', 'AC---GT', 'ACAAAGT'),
            ('ACAAAAGT', 'AC----GT', ''),
            ('ACAAACCCGT', 'AC------GT', 'ACAAACCCGT'),

            # both
            ('AG-GTTT', 'ACC-TTT', 'AGXTTT'),
        )
        for seq, ref, expected in cases:
            self._test(seq, ref, expected)

    def test_mismatched_lengths(self):
        seq = 'ACC'
        ref = 'AC'
        self.assertRaises(ValueError, correct_shifts, seq, ref)


unittest.main()
