import pickle
import unittest
import sys

import pyopal


class TestScoreMatrix(unittest.TestCase):

    def test_aa_blosum50(self):
        aa = pyopal.ScoreMatrix.aa("BLOSUM50")
        matrix = aa.matrix
        diagonal = [ matrix[i][i] for i in range(len(aa.alphabet)) ]
        self.assertEqual(diagonal, [5, 7, 7, 8, 13, 7, 6, 8, 10, 5, 5, 6, 7, 8, 10, 5, 5, 15, 8, 5, 6, 5 ,-1 ,1])

    def test_aa_blosum62(self):
        aa = pyopal.ScoreMatrix.aa("BLOSUM62")
        matrix = aa.matrix
        diagonal = [ matrix[i][i] for i in range(len(aa.alphabet)) ]
        self.assertEqual(diagonal, [4, 5, 6, 6, 9, 5, 5, 6, 8, 4, 4, 5, 5, 6, 7, 4, 5, 11, 7, 4, 4, 4, -1])

    def test_aa_invalid_name(self):
        with self.assertRaises(ValueError):
            aa = pyopal.ScoreMatrix.aa("nonsensical")

    @unittest.skipUnless(sys.implementation.name == "cpython", "memoryview not supported")
    def test_memoryview(self):
        aa = pyopal.ScoreMatrix.aa("BLOSUM50")
        mem = memoryview(aa)

        self.assertEqual(mem.shape, (24, 24))
        self.assertEqual(mem[0, 0], 5)
        self.assertEqual(mem[3, 3], 8)
        
        l = mem.tolist()
        self.assertListEqual(l, aa.matrix)

    def test_init_invalid_length(self):
        with self.assertRaises(ValueError):
            aa = pyopal.ScoreMatrix(
                "ATGC",
                [
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                ]
            )
        with self.assertRaises(ValueError):
            aa = pyopal.ScoreMatrix(
                "ATGC",
                [
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0],
                ]
            )

    def test_eq(self):
        sm1 = pyopal.ScoreMatrix.aa("BLOSUM50")
        sm2 = pyopal.ScoreMatrix.aa("BLOSUM50")
        sm3 = pyopal.ScoreMatrix.aa("BLOSUM62")
        self.assertEqual(sm1, sm1)
        self.assertEqual(sm1, sm2)
        self.assertNotEqual(sm1, sm3)
        self.assertNotEqual(sm1, 12)

    def test_pickle(self):
        sm1 = pyopal.ScoreMatrix.aa()
        sm2 = pickle.loads(pickle.dumps(sm1))
        self.assertEqual(sm1.alphabet, sm2.alphabet)
        self.assertEqual(sm1.matrix, sm2.matrix)