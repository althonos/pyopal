import unittest

import pyopal


class TestSearch(unittest.TestCase):
    def test_test1_sw(self):
        # #0: 47 (0, 1) (5, 7)
        # T: ACCGCTG (1 - 7)
        # Q: ACCTC_G (0 - 5)

        query = "ACCTCG"
        target = "AACCGCTG"
        db = pyopal.Database([target])

        results = db.search(query, algorithm="sw")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.ScoreResult)
        self.assertEqual(results[0].score, 47)

        results = db.search(query, algorithm="sw", mode="score")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.ScoreResult)
        self.assertEqual(results[0].score, 47)

        results = db.search(query, algorithm="sw", mode="end")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.EndResult)
        self.assertEqual(results[0].score, 47)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].target_end, 7)

        results = db.search(query, algorithm="sw", mode="full")
        self.assertEqual(len(db), 1)
        self.assertEqual(results[0].score, 47)
        self.assertIsInstance(results[0], pyopal.FullResult)
        self.assertIsNot(results[0].alignment, None)
        self.assertEqual(results[0].query_start, 0)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].target_start, 1)
        self.assertEqual(results[0].target_end, 7)

        self.assertAlmostEqual(results[0].coverage("query"), 1)
        self.assertAlmostEqual(results[0].coverage("target"), 7/8)

    def test_test1_nw(self):
        # #0: 44 (0, 0) (5, 7)
        # T: AACCGCTG (0 - 7)
        # Q: _ACCTC_G (0 - 5)

        query = "ACCTCG"
        target = "AACCGCTG"
        db = pyopal.Database([target])

        results = db.search(query, algorithm="nw")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.ScoreResult)
        self.assertEqual(results[0].score, 44)

        results = db.search(query, algorithm="nw", mode="score")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.ScoreResult)
        self.assertEqual(results[0].score, 44)

        results = db.search(query, algorithm="nw", mode="end")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.EndResult)
        self.assertEqual(results[0].score, 44)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].target_end, 7)

        results = db.search(query, algorithm="nw", mode="full")
        self.assertEqual(len(db), 1)
        self.assertEqual(results[0].score, 44)
        self.assertIsInstance(results[0], pyopal.FullResult)
        self.assertIsNot(results[0].alignment, None)
        self.assertEqual(results[0].query_start, 0)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].target_start, 0)
        self.assertEqual(results[0].target_end, 7)
        self.assertEqual(results[0].coverage("query"), 1)
        self.assertEqual(results[0].coverage("target"), 7/8)