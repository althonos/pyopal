# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True
"""Bindings to Opal, a SIMD-accelerated pairwise sequence aligner.

References:
    - Korpar M., Šošić M., Blažeka D., Šikić M.
      *SW#db: GPU-Accelerated Exact Sequence Similarity Database Search*.
      PLoS One. 2015;10(12):e0145857. Published 2015 Dec 31.
      :doi:`10.1371/journal.pone.0145857`.

"""

# --- C imports ----------------------------------------------------------------

from libc.limits cimport UCHAR_MAX
from libcpp cimport bool, nullptr
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr

from scoring_matrices cimport ScoringMatrix

from . cimport opal
from .opal cimport OpalSearchResult
from .shared_mutex cimport shared_mutex

# --- Constants ----------------------------------------------------------------

cdef extern from * nogil:
    """
    const size_t MAX_ALPHABET_SIZE = 32;
    """
    const size_t MAX_ALPHABET_SIZE

# --- Type definitions ---------------------------------------------------------

ctypedef unsigned char       digit_t
ctypedef shared_ptr[digit_t] seq_t
ctypedef OpalSearchResult*   OpalSearchResult_ptr

ctypedef int (*searchfn_t)(
    unsigned char*,
    int,
    unsigned char**,
    int,
    int*,
    int,
    int,
    int*,
    int,
    OpalSearchResult**,
    const int,
    int,
    int,
) noexcept nogil

# --- Exclusive read/write lock ------------------------------------------------

cdef class SharedMutex:
    cdef          shared_mutex mutex
    cdef readonly ReadLock     read
    cdef readonly WriteLock    write

cdef class ReadLock:
    cdef SharedMutex owner

cdef class WriteLock:
    cdef SharedMutex owner

# --- Python classes -----------------------------------------------------------

cdef class Alphabet:
    cdef readonly str   letters
    cdef readonly int   length
    cdef char          _unknown
    cdef unsigned char _letters[MAX_ALPHABET_SIZE]
    cdef char          _ahash[UCHAR_MAX]

    cpdef void encode_into(self, const unsigned char[:] sequence, digit_t[:] encoded)
    cpdef void decode_into(self, const digit_t[:] encoded, unsigned char[:] sequence)

    cpdef bytes encode(self, object sequence)
    cpdef str decode(self, object encoded)

# --- Sequence storage ---------------------------------------------------------

cdef class BaseDatabase:
    cdef readonly SharedMutex      lock
    cdef readonly Alphabet         alphabet

    cdef const digit_t** get_sequences(self) except? NULL
    cdef const int*      get_lengths(self) except? NULL
    cdef size_t          get_size(self) noexcept


cdef class Database(BaseDatabase):
    cdef vector[seq_t]    _sequences
    cdef vector[digit_t*] _pointers
    cdef vector[int]      _lengths

    cdef digit_t* _encode(self, object sequence) except *

    cpdef void clear(self) except *
    cpdef void extend(self, object sequences) except *
    cpdef void append(self, object sequence) except *
    cpdef void reverse(self) except *
    cpdef void insert(self, ssize_t index, object sequence)

    cpdef Database mask(self, object bitmask)
    cpdef Database extract(self, object indices)

# --- Aligner ------------------------------------------------------------------

cdef class ScoreResult:
    cdef ssize_t          _target_index
    cdef OpalSearchResult _result

cdef class EndResult(ScoreResult):
    pass

cdef class FullResult(EndResult):
    cdef int _query_length
    cdef int _target_length

    cpdef str cigar(self)
    cpdef float identity(self)
    cpdef float coverage(self, str reference=*)

cdef class Aligner:
    cdef readonly Alphabet      alphabet
    cdef readonly ScoringMatrix scoring_matrix
    cdef readonly int           gap_open
    cdef readonly int           gap_extend
    
    cdef          searchfn_t    _search
    cdef          int*          _int_matrix