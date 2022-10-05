import unittest

import pyopal


class TestSearch(unittest.TestCase):
    def test_test1(self):

        query = "ACCTCG"
        target = "AACCGCTG"

        db = pyopal.Database([target])

        results = db.search(query)
        self.assertEqual(len(db), 1)
        self.assertEqual(results[0].score, 47)
        self.assertIs(results[0].alignment, None)
        self.assertIs(results[0].query_start, None)
        self.assertIs(results[0].query_end, None)
        self.assertIs(results[0].target_start, None)
        self.assertIs(results[0].target_end, None)

        results = db.search(query)
        self.assertEqual(len(db), 1)
        self.assertEqual(results[0].score, 47)
        self.assertIs(results[0].alignment, None)
        self.assertIs(results[0].query_start, None)
        self.assertIs(results[0].query_end, None)
        self.assertIs(results[0].target_start, None)
        self.assertIs(results[0].target_end, None)

        results = db.search(query, mode="end")
        self.assertEqual(len(db), 1)
        self.assertEqual(results[0].score, 47)
        self.assertIs(results[0].alignment, None)
        self.assertIs(results[0].query_start, None)
        self.assertEqual(results[0].query_end, 5)
        self.assertIs(results[0].target_start, None)
        self.assertEqual(results[0].target_end, 7)

        results = db.search(query, mode="full")
        self.assertEqual(len(db), 1)
        self.assertEqual(results[0].score, 47)
        self.assertIsNot(results[0].alignment, None)
        self.assertEqual(results[0].query_start, 0)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].target_start, 1)
        self.assertEqual(results[0].target_end, 7)
