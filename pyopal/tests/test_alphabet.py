import pickle
import unittest

import pyopal


class TestAlphabet(unittest.TestCase):

    def test_len(self):
        alphabet = pyopal.Alphabet()
        self.assertEqual(len(alphabet), 24)
        self.assertEqual(len(alphabet), len(alphabet.letters))

        alphabet = pyopal.Alphabet("ATGC")
        self.assertEqual(len(alphabet), 4)
        self.assertEqual(len(alphabet), len(alphabet.letters))

    def test_letters(self):
        alphabet = pyopal.Alphabet("ATGC")
        self.assertEqual(alphabet.letters, "ATGC")

    def test_init_error_duplicate_letters(self):
        with self.assertRaises(ValueError):
            alphabet = pyopal.Alphabet("AAAA")

    def test_init_error_lowercase_letters(self):
        with self.assertRaises(ValueError):
            alphabet = pyopal.Alphabet("AtgC")

    def test_init_error_invalid_symbols(self):
        with self.assertRaises(ValueError):
            alphabet = pyopal.Alphabet("A[]C")

    def test_encode_str(self):
        alphabet = pyopal.Alphabet("ATGC")

        s1 = alphabet.encode("ATGC")
        self.assertEqual(s1, bytes([0, 1, 2, 3]))

        s2 = alphabet.encode("AAAAA")
        self.assertEqual(s2, bytes([0, 0, 0, 0, 0]))

    def test_encode_bytes(self):
        alphabet = pyopal.Alphabet("ATGC")

        s1 = alphabet.encode(b"ATGC")
        self.assertEqual(s1, bytes([0, 1, 2, 3]))

        s2 = alphabet.encode(b"AAAAA")
        self.assertEqual(s2, bytes([0, 0, 0, 0, 0]))

    def test_decode_bytes(self):
        alphabet = pyopal.Alphabet("ATGC")

        s1 = alphabet.decode(bytes([0, 1, 2, 3]))
        self.assertEqual(s1, "ATGC")

        s1 = alphabet.decode(bytes([0, 0, 0, 0, 0]))
        self.assertEqual(s1, "AAAAA")

    def test_decode_bytearray(self):
        alphabet = pyopal.Alphabet("ATGC")

        s1 = alphabet.decode(bytearray([0, 1, 2, 3]))
        self.assertEqual(s1, "ATGC")

        s1 = alphabet.decode(bytearray([0, 0, 0, 0, 0]))
        self.assertEqual(s1, "AAAAA")

    def test_decode_memoryview(self):
        alphabet = pyopal.Alphabet("ATGC")

        s1 = alphabet.decode(memoryview(bytearray([0, 1, 2, 3])))
        self.assertEqual(s1, "ATGC")

        s1 = alphabet.decode(memoryview(bytearray([0, 0, 0, 0, 0])))
        self.assertEqual(s1, "AAAAA")