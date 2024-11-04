import pickle
import unittest
import sys

import pyopal


class TestScoreResult(unittest.TestCase):

    def test_init(self):
        result = pyopal.ScoreResult(10, score=30)
        self.assertEqual(result.score, 30)
        self.assertEqual(result.target_index, 10)

    def test_repr(self):
        result = pyopal.ScoreResult(target_index=10, score=30)
        self.assertEqual(repr(result), "ScoreResult(10, score=30)")

    def test_pickle(self):
        r1 = pyopal.ScoreResult(target_index=10, score=30)
        r2 = pickle.loads(pickle.dumps(r1))
        self.assertEqual(r2.score, 30)
        self.assertEqual(r2.target_index, 10)

    def test_eq(self):
        r1 = pyopal.ScoreResult(target_index=10, score=30)
        r2 = pyopal.ScoreResult(target_index=10, score=30)
        r3 = pyopal.ScoreResult(target_index=12, score=50)
        self.assertEqual(r1, r1)
        self.assertEqual(r1, r2)
        self.assertNotEqual(r1, r3)
        self.assertNotEqual(r1, 12)


class TestEndResult(unittest.TestCase):

    def test_init(self):
        result = pyopal.EndResult(2, score=30, query_end=10, target_end=20)
        self.assertEqual(result.score, 30)
        self.assertEqual(result.target_index, 2)
        self.assertEqual(result.query_end, 10)
        self.assertEqual(result.target_end, 20)

    def test_repr(self):
        result = pyopal.EndResult(target_index=10, score=30, query_end=10, target_end=20)
        self.assertEqual(repr(result), "EndResult(10, score=30, query_end=10, target_end=20)")

    def test_pickle(self):
        r1 = pyopal.EndResult(target_index=10, score=30, query_end=10, target_end=20)
        r2 = pickle.loads(pickle.dumps(r1))
        self.assertEqual(r2.score, 30)
        self.assertEqual(r2.target_index, 10)
        self.assertEqual(r2.query_end, 10)
        self.assertEqual(r2.target_end, 20)

    def test_eq(self):
        r1 = pyopal.EndResult(target_index=10, score=30, query_end=10, target_end=20)
        r2 = pyopal.EndResult(target_index=10, score=30, query_end=10, target_end=20)
        r3 = pyopal.EndResult(target_index=10, score=35, query_end=20, target_end=60)
        self.assertEqual(r1, r1)
        self.assertEqual(r1, r2)
        self.assertNotEqual(r1, r3)
        self.assertNotEqual(r1, 12)
        

class TestFullResult(unittest.TestCase):

    def test_init(self):
        result = pyopal.FullResult(10, score=30, query_end=10, target_end=20, query_start=0, target_start=10, query_length=100, target_length=100, alignment="M"*10)
        self.assertEqual(result.score, 30)
        self.assertEqual(result.target_index, 10)
        self.assertEqual(result.query_end, 10)
        self.assertEqual(result.target_end, 20)
        self.assertEqual(result.query_start, 0)
        self.assertEqual(result.target_start, 10)
        self.assertEqual(result.query_length, 100)
        self.assertEqual(result.target_length, 100)
        self.assertEqual(result.alignment, "M"*10)

    def test_pickle(self):
        r1 = pyopal.FullResult(10, score=30, query_end=10, target_end=20, query_start=0, target_start=10, query_length=100, target_length=100, alignment="M"*10)
        r2 = pickle.loads(pickle.dumps(r1))
        self.assertEqual(r2.score, 30)
        self.assertEqual(r2.target_index, 10)
        self.assertEqual(r2.query_end, 10)
        self.assertEqual(r2.target_end, 20)
        self.assertEqual(r2.query_start, 0)
        self.assertEqual(r2.target_start, 10)
        self.assertEqual(r2.query_length, 100)
        self.assertEqual(r2.target_length, 100)
        self.assertEqual(r2.alignment, "M"*10)

    def test_eq(self):
        r1 = pyopal.FullResult(10, score=30, query_end=10, target_end=20, query_start=0, target_start=10, query_length=100, target_length=100, alignment="M"*10)
        r2 = pyopal.FullResult(10, score=30, query_end=10, target_end=20, query_start=0, target_start=10, query_length=100, target_length=100, alignment="M"*10)
        r3 = pyopal.FullResult(2, score=48, query_end=10, target_end=20, query_start=0, target_start=30, query_length=500, target_length=200, alignment="M"*10)
        self.assertEqual(r1, r1)
        self.assertEqual(r1, r2)
        self.assertNotEqual(r1, r3)
        self.assertNotEqual(r1, 12)