import unittest

import pyopal


class TestDatabase(unittest.TestCase):
    def test_getitem(self):
        sequences = ["ATGC", "ATTC", "TTCG"]
        database = pyopal.Database(sequences)
        self.assertEqual(database[0], sequences[0])
        self.assertEqual(database[1], sequences[1])
        self.assertEqual(database[2], sequences[2])
        self.assertEqual(database[-1], sequences[-1])
        self.assertEqual(database[-2], sequences[-2])
        self.assertEqual(database[-3], sequences[-3])

        sequences = ["ATGC", "ATTC", "TTCG"]
        database = pyopal.Database([seq.encode("ascii") for seq in sequences])
        self.assertEqual(database[0], sequences[0])
        self.assertEqual(database[1], sequences[1])
        self.assertEqual(database[2], sequences[2])
        self.assertEqual(database[-1], sequences[-1])
        self.assertEqual(database[-2], sequences[-2])
        self.assertEqual(database[-3], sequences[-3])

    def test_getitem_index_error(self):
        sequences = ["ATGC", "ATTC", "TTCG"]
        database = pyopal.Database(sequences)
        with self.assertRaises(IndexError):
            database[3]
        with self.assertRaises(IndexError):
            database[-4]
        with self.assertRaises(IndexError):
            database[-8]

    def test_reverse(self):
        sequences = ["ATGC", "ATTC", "TTCG"]
        database = pyopal.Database(sequences)
        self.assertEqual(list(database), sequences)
        database.reverse()
        self.assertEqual(list(database), list(reversed(sequences)))

    def test_reverse_empty(self):
        database = pyopal.Database()
        self.assertEqual(len(database), 0)
        database.reverse()
        self.assertEqual(len(database), 0)
