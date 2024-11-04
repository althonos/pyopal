import random
import unittest

import pyopal
import multiprocessing.pool

class TestAlignThreads(unittest.TestCase):

    def test_test1(self):
        # #0: 44 (0, 0) (5, 7)
        # T: AACCGCTG (0 - 7)
        # Q: _ACCTC_G (0 - 5)

        query = "ACCTCG"
        target = ["AACCGCTG", "AACCGCTA", "AACCGCTC", "AACCGCTT"]
        results = list(pyopal.align(query, target, threads=1, mode="full", algorithm="nw", ordered=True))
        self.assertEqual(results[0].target_index, 0)
        self.assertEqual(results[0].target_start, 0)
        self.assertEqual(results[0].target_end, 7)
        self.assertEqual(results[0].query_start, 0)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].score, 44)

    def test_2(self):
        # #0: 44 (0, 0) (5, 7)
        # T: AACCGCTG (0 - 7)
        # Q: _ACCTC_G (0 - 5)

        query = "ACCTCG"
        target = ["AACCGCTG", "AACCGCTA", "AACCGCTC", "AACCGCTT"]
        results = list(pyopal.align(query, target, threads=2, mode="full", algorithm="nw", ordered=True))
        self.assertEqual(results[0].target_index, 0)
        self.assertEqual(results[0].target_start, 0)
        self.assertEqual(results[0].target_end, 7)
        self.assertEqual(results[0].query_start, 0)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].score, 44)