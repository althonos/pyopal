from libcpp.vector cimport vector

cdef extern from "ScoreMatrix.hpp" nogil:

    cdef cppclass ScoreMatrix:
        ScoreMatrix()
        ScoreMatrix(vector[unsigned char] alphabet, vector[int] matrix)
        ScoreMatrix(const char* filepath) except +

        int getAlphabetLength()
        unsigned char* getAlphabet()
        int* getMatrix()

        @staticmethod
        ScoreMatrix getBlosum50()