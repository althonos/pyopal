from libc.stdint cimport uint32_t
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr

from sword.chain cimport Chain, ChainSet

cdef extern from "reader.hpp" nogil:

    cppclass Reader:
        bool read_chains(ChainSet& dst, size_t max_bytes)
        unique_ptr[Reader] createReader(const string& path)

    unique_ptr[Reader] createReader(const string& path)