from libc.stdint cimport uint32_t
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr

from sword.score_matrix cimport ScoreMatrix

cdef extern from "kmers.hpp" nogil:

    cppclass Kmers:
        uint32_t kmer_length() const

        const vector[uint32_t]& kmer_substitutions(uint32_t kmer) const
        const vector[uint32_t]& kmer_substitutions(const string& kmer) const

    unique_ptr[Kmers] createKmers(uint32_t kmer_length, uint32_t score_threshold, shared_ptr[ScoreMatrix] score_matrix)