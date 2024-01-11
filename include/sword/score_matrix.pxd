from libc.stdint cimport int32_t, uint32_t
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr

cdef extern from "score_matrix.hpp" nogil:

    cppclass ScoreMatrixType:
        pass

cdef extern from "score_matrix.hpp" namespace "ScoreMatrixType" nogil:
    const ScoreMatrixType kBlosum45
    const ScoreMatrixType kBlosum50
    const ScoreMatrixType kBlosum62
    const ScoreMatrixType kBlosum80
    const ScoreMatrixType kBlosum90
    const ScoreMatrixType kPam30
    const ScoreMatrixType kPam70
    const ScoreMatrixType kPam250


cdef extern from "score_matrix.hpp" nogil:

    cppclass ScoreMatrix:
        int* data()
        const int* data() const
        ScoreMatrixType type() const
        int32_t gap_open() const
        int32_t gap_extend() const

        string scorerName()

        int score(uint32_t row, uint32_t column) const

    unique_ptr[ScoreMatrix] createScoreMatrix(ScoreMatrixType type, int32_t gap_open, int32_t gap_extend)


cdef extern from "score_matrix.hpp" namespace "ScoreMatrix" nogil:
    const uint32_t num_rows_
    const uint32_t num_columns_
    const uint32_t size_