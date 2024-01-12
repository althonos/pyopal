from libc.stdint cimport uint32_t
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr

cdef extern from "reader.hpp" nogil:
    cppclass Reader:
        pass

cdef extern from "chain.hpp" nogil:

    cppclass Chain:
        const uint32_t id() const
        const string& name() const
        const size_t name_length() const
        const string& data() const
        const size_t length() const

        unique_ptr[Chain] createChain(uint32_t id, char* name, uint32_t name_length, char* data, uint32_t data_length)

    ctypedef vector[unique_ptr[Chain]] ChainSet 
    unique_ptr[Chain] createChain(uint32_t id, char* name, uint32_t name_length, char* data, uint32_t data_length)
    void createChainSet(ChainSet& dst, const string& path)

    unique_ptr[Reader] createChainSetPartInitialize(string& path)
    bool createChainSetPart(ChainSet& dst, shared_ptr[Reader] reader, size_t max_bytes)