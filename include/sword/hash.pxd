from libc.stdint cimport uint32_t
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr, unique_ptr

from sword.kmers cimport Kmers
from sword.chain cimport ChainSet


cdef extern from "hash.hpp" nogil:

    cppclass Hit:
        Hit()
        Hit(uint32_t id, uint32_t position)

        uint32_t id() const
        uint32_t position() const


cdef extern from "hash.hpp" namespace "Hash" nogil:

    ctypedef vector[Hit].iterator Iterator


cdef extern from "hash.hpp" nogil:

    cppclass Hash:
        void hits(Iterator& start, Iterator& end, uint32_t key)

    unique_ptr[Hash] createHash(const ChainSet& chains, uint32_t start, uint32_t length, shared_ptr[Kmers] kmers)