# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True

cimport libcpp.algorithm
from libc.stdint cimport int32_t, uint32_t, uint64_t
from libcpp cimport bool, nullptr
from libcpp.memory cimport unique_ptr, shared_ptr

cimport sword.kmers
cimport sword.reader
cimport sword.score_matrix
from sword.kmers cimport Kmers as _Kmers
from sword.chain cimport ChainSet as _ChainSet, Chain as _Chain
from sword.reader cimport Reader as _Reader 
from sword.score_matrix cimport ScoreMatrix as _ScoreMatrix, ScoreMatrixType as _ScoreMatrixType

import os


cdef class Kmers:
    cdef shared_ptr[_Kmers] _kmers

    def __init__(self, ScoreMatrix score_matrix, kmer_length = 3, score_threshold = 13):
        self._kmers = shared_ptr[_Kmers](
            sword.kmers.createKmers(kmer_length, score_threshold, score_matrix._score_matrix)
        )


cdef class ChainSet:
    cdef _ChainSet _chains

    def __len__(self):
        return self._chains.size()

    def __getitem__(self, ssize_t i):
        cdef ssize_t _i = i
        if _i < 0:
            _i += self._chains.size()
        if _i < 0 or _i >= self._chains.size():
            raise IndexError(i)
        return self._chains[i].get().name()


cdef class Reader:
    cdef unique_ptr[_Reader] reader

    def __init__(self, path):
        cdef bytes _path = os.fsencode(path)
        self.reader = sword.reader.createReader(_path)

    cpdef ChainSet read(self):
        assert self._reader
        cdef ChainSet s = ChainSet()
        self.reader.get().read_chains(s._chains, 0)
        return s


cdef class ScoreMatrix:
    cdef shared_ptr[_ScoreMatrix] _score_matrix

    def __init__(self, str name = "BLOSUM62", int32_t gap_open = 10, int32_t gap_extend = 1):
        cdef _ScoreMatrixType ty
        if name == "BLOSUM62":
            ty = sword.score_matrix.kBlosum62
        else:
            raise ValueError(f"unsupported score matrix: {name!r}")
        self._score_matrix = shared_ptr[_ScoreMatrix](
            sword.score_matrix.createScoreMatrix(
                ty,
                gap_open,
                gap_extend,
            )
        )

    @property
    def name(self):
        assert self._score_matrix != nullptr
        return self._score_matrix.get().scorerName().decode()


cdef extern from *:
    """
    bool chainLengthKey(const std::unique_ptr<Chain>& left, const std::unique_ptr<Chain>& right) {
        return left->length() < right->length();
    }
    """
    bool chainLengthKey(const unique_ptr[_Chain]& left, const unique_ptr[_Chain]& right)


cpdef object preprocess_database(
    ChainSet database,
    size_t num_threads = 1,
    uint32_t max_short_length = 2000,
):
    # ported from: preprocDatabase in `database_search.cpp`
    cdef vector[uint32_t] dst

    # sort by length
    libcpp.algorithm.sort(database._chains.begin(), database._chains.end(), chainLengthKey)

    # split tasks between long and short
    cdef uint64_t short_total_length = 0
    cdef uint64_t long_total_length  = 0
    cdef uint32_t split              = 0

    for i in range(len(database)):
        l = database._chains[i].length() 
        if l > max_short_length:
            if split == 0:
                split = i
            long_total_length += l
        else:
            short_total_length += l

    if short_total_length == 0:
        split = 0
    if long_total_length == 0:
        split = len(database)


    # spread tasks across threads
    cdef uint64_t short_task_size = short_total_length / num_threads
    cdef uint64_t long_task_size = long_total_length / num_threads

    dst.reserve(2*num_threads + 1)

    #
    cdef uint64_t total_length = 0



cpdef object search_database(
    str database_path,
    str queries_path,
    uint32_t kmer_length = 3,
    uint32_t max_candidates = 30000,
    ScoreMatrix score_matrix = ScoreMatrix(),
    uint32_t score_threshold = 13,
):
    cdef ChainSet queries
    cdef ChainSet db      = ChainSet()
    cdef Reader   dreader = Reader(database_path)
    cdef Reader   qreader = Reader(queries_path)

    # read queries
    queries = qreader.read()
    libcpp.algorithm.sort(queries._chains.begin(), queries._chains.end(), chainLengthKey)

    # create initial kmers
    cdef Kmers kmers = Kmers(score_matrix, kmer_length, score_threshold)

    while True:

        sword.chain.createChainSetPart(db._chain, dreader._reader)







