import random
import unittest

import pyopal
from scoring_matrices import ScoringMatrix


class TestAligner(unittest.TestCase):

    def test_init_scoring_matrix_str(self):
        matrix = ScoringMatrix.from_name("VTML80")
        aligner = pyopal.Aligner("VTML80")
        self.assertEqual(aligner.scoring_matrix, matrix)

    def test_init_scoring_matrix_object(self):
        matrix = ScoringMatrix.from_name("VTML80")
        aligner = pyopal.Aligner(matrix)
        self.assertEqual(aligner.scoring_matrix, matrix)

    def test_init_scoring_matrix_error(self):
        self.assertRaises(TypeError, pyopal.Aligner, 1)


class _TestSearchOverflow(object):

    algorithm = NotImplemented

    def test_overflow(self):
        r = random.Random(0)
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        proteins = [
            "".join(random.choices(alphabet, k=k))
            for k in range(1000, 36000, 1000)
        ]
        database = pyopal.Database(proteins)
        aligner = pyopal.Aligner()
        results = aligner.align(proteins[0], database, mode="score", algorithm=self.algorithm)

class TestSearchNW(unittest.TestCase, _TestSearchOverflow):
    algorithm = "nw"

    def test_test1(self):
        # #0: 44 (0, 0) (5, 7)
        # T: AACCGCTG (0 - 7)
        # Q: _ACCTC_G (0 - 5)

        query = "ACCTCG"
        target = "AACCGCTG"
        aligner = pyopal.Aligner()
        db = pyopal.Database([target])

        results = aligner.align(query, db, algorithm=self.algorithm)
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.ScoreResult)
        self.assertEqual(results[0].score, 44)

        results = aligner.align(query, db, algorithm=self.algorithm, mode="score")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.ScoreResult)
        self.assertEqual(results[0].score, 44)

        results = aligner.align(query, db, algorithm=self.algorithm, mode="end")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.EndResult)
        self.assertEqual(results[0].score, 44)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].target_end, 7)

        results = aligner.align(query, db, algorithm=self.algorithm, mode="full")
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


class TestSearchOV(unittest.TestCase, _TestSearchOverflow):
    algorithm = "ov"


class TestSearchHW(unittest.TestCase, _TestSearchOverflow):
    algorithm = "hw"


class TestSearchSW(unittest.TestCase, _TestSearchOverflow):
    algorithm = "sw"

    def test_test1(self):
        # #0: 47 (0, 1) (5, 7)
        # T: ACCGCTG (1 - 7)
        # Q: ACCTC_G (0 - 5)

        query = "ACCTCG"
        target = "AACCGCTG"
        aligner = pyopal.Aligner()
        db = pyopal.Database([target])

        results = aligner.align(query, db, algorithm=self.algorithm)
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.ScoreResult)
        self.assertEqual(results[0].score, 47)

        results = aligner.align(query, db, algorithm=self.algorithm, mode="score")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.ScoreResult)
        self.assertEqual(results[0].score, 47)

        results = aligner.align(query, db, algorithm=self.algorithm, mode="end")
        self.assertEqual(len(db), 1)
        self.assertIsInstance(results[0], pyopal.EndResult)
        self.assertEqual(results[0].score, 47)
        self.assertEqual(results[0].query_end, 5)
        self.assertEqual(results[0].target_end, 7)

        results = aligner.align(query, db, algorithm=self.algorithm, mode="full")
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

