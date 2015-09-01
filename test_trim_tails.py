import unittest
import numpy as np
from trim_tails import trim_tail


class TestTrimTails(unittest.TestCase):

    def test_null_case(self):
        seq = "ACGT" * 5
        result, _ = trim_tail(seq, "A", 0.25, 0.99)
        self.assertEquals(result, seq)

    def test_simple_case(self):
        head = "ACGT" * 5
        tail = "A" * 10
        seq = "".join([head, tail])
        result, _ = trim_tail(seq, "A", 0.25, 0.99)
        self.assertEquals(result, head)

    def test_rev_case(self):
        head = "T" * 10
        tail = "ACGT" * 5
        seq = "".join([head, tail])
        result, _ = trim_tail(seq, "T", 0.25, 0.99, reverse=True)
        self.assertEquals(result, tail)

    def test_with_error(self):
        head = "ACGT" * 5
        tail = "A" * 5 + "G" + "A" * 5
        seq = "".join([head, tail])
        result, _ = trim_tail(seq, "A", 0.25, 0.99)
        self.assertEquals(result, head)

    def test_penalty(self):
        head = "ACGT" * 5
        tail = "A" * 10
        seq = "".join([head, tail])
        result, _ = trim_tail(seq, "A", 0.25, 0.99, penalty=np.log(0.000001))
        self.assertEquals(result, seq)
        

if __name__ == '__main__':
    unittest.main()
