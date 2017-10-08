import unittest

from flea.correct_shifts import correct_shifts

class TestCorrectShifts(unittest.TestCase):

    def _test(self, seq, ref, expected, keep=False, deletion_strategy='discard'):
        self.assertEquals(correct_shifts(seq, ref, keep=keep,
                                         deletion_strategy=deletion_strategy)[0], expected)

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

            # both
            ('AG-GTTT', 'ACC-TTT', ''),
        )
        for seq, ref, expected in cases:
            self._test(seq, ref, expected)

    def test_keep(self):
        cases = (
            ('ACAAAAG', 'AC----G', 'ACAAAAG'),

            # deletions
            ('AC-ACG', 'ACGACG', 'ACNACG'),
            ('AC--CG', 'ACGACG', 'ACNNCG'),
            ('AC---T', 'ACGACG', 'ACT'),
            ('AC----ACG', 'ACGACGACG', 'ACNACG'),

            # both; keep=True
            ('AC-TTTT', 'ACCTTT-', 'ACNTTT'),
        )
        for seq, ref, expected in cases:
            self._test(seq, ref, expected, keep=True,
                       deletion_strategy='n')

    def test_correct(self):
        cases = (
            # deletions
            ('AC-ACG', 'ACGACG', 'ACGACG'),
        )
        for seq, ref, expected in cases:
            self._test(seq, ref, expected, deletion_strategy='reference')

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
