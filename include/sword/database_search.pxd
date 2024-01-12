from libc.stdint cimport uint32_t
from libcpp.vector cimport vector

cdef extern from * nogil:
    """
    class ChainEntry {
    public:
        ChainEntry() {};
        ChainEntry(uint32_t chain_idx, uint32_t data)
                : chain_idx_(chain_idx), data_(data) {
        };

        ~ChainEntry() = default;

        uint32_t chain_idx() const {
            return chain_idx_;
        }

        uint32_t data() const {
            return data_;
        }

    private:

        uint32_t chain_idx_;
        uint32_t data_;
    };

    using Indexes = std::vector<std::vector<uint32_t>>;
    using ChainEntrySet = std::vector<std::vector<ChainEntry>>;
    """
    cppclass ChainEntry:
        ChainEntry()
        ChainEntry(uint32_t chain_idx, uint32_t data)

        uint32_t chain_idx() const
        uint32_t data() const

    ctypedef vector[vector[ChainEntry]] ChainEntrySet

#cdef extern from "database_search.hpp" nogil:
    ctypedef vector[vector[uint32_t]] Indexes

