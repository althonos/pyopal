import pickle
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

    def test_pickle(self):
        sequences = ["ATGC", "ATTC", "TTCG"]
        database = pyopal.Database(sequences)
        unpickled = pickle.loads(pickle.dumps(database))
        self.assertEqual(list(unpickled), sequences)

    def test_insert(self):
        sequences = ["ATGC", "ATTC"]
        database = pyopal.Database(sequences)
        database.insert(1, "TTCC")
        self.assertEqual(list(database), ["ATGC", "TTCC", "ATTC"])
        database.insert(-10, "TTTT")
        self.assertEqual(list(database), ["TTTT", "ATGC", "TTCC", "ATTC"])
        database.insert(10, "AAAA")
        self.assertEqual(list(database), ["TTTT", "ATGC", "TTCC", "ATTC", "AAAA"])

    def test_delitem(self):
        sequences = ["ATGC", "ATTC", "TTCG"]
        database = pyopal.Database(sequences)
        self.assertEqual(list(database), sequences)
        del database[1]
        self.assertEqual(list(database), ["ATGC", "TTCG"])
        del database[-2]
        self.assertEqual(list(database), ["TTCG"])
        del database[0]
        self.assertEqual(list(database), [])
        with self.assertRaises(IndexError):
            del database[0]
        with self.assertRaises(IndexError):
            del database[-1]

    def test_setitem(self):
        sequences = ["ATGC", "ATTC", "TTCG"]
        database = pyopal.Database(sequences)
        self.assertEqual(list(database), sequences)
        database[2] = "AAAT"
        self.assertEqual(list(database), ["ATGC", "ATTC", "AAAT"])
        with self.assertRaises(IndexError):
            database[-8] = "TCGA"
        with self.assertRaises(IndexError):
            database[5] = "TCGA"
