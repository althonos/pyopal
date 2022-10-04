# distutils: language = c++
# cython: language_level=3, linetrace=True, embedsignature=True, binding=True
"""Bindings to Opal, a SIMD-accelerated pairwise sequence aligner.

References:
    - Korpar M., Šošić M., Blažeka D., Šikić M. 
      *SW#db: GPU-Accelerated Exact Sequence Similarity Database Search*. 
      PLoS One. 2015;10(12):e0145857. Published 2015 Dec 31. 
      :doi:`10.1371/journal.pone.0145857`.

"""

# --- C imports ----------------------------------------------------------------

from cpython.bytes cimport PyBytes_AsString, PyBytes_FromStringAndSize
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from libc.string cimport memset

cimport opal
from opal cimport OpalSearchResult
from opal.score_matrix cimport ScoreMatrix

cdef extern from "<cctype>" namespace "std" nogil:
    cdef int toupper(int ch)
    cdef int tolower(int ch)


# --- Runtime CPU detection ----------------------------------------------------

_TARGET_CPU           = TARGET_CPU
_SSE2_RUNTIME_SUPPORT = False
_SSE2_BUILD_SUPPORT   = False
_SSE4_RUNTIME_SUPPORT = False
_SSE4_BUILD_SUPPORT   = False
_AVX2_RUNTIME_SUPPORT = False
_AVX2_BUILD_SUPPORT   = False
_NEON_RUNTIME_SUPPORT = False
_NEON_BUILD_SUPPORT   = False

IF TARGET_CPU == "x86" and TARGET_SYSTEM in ("freebsd", "linux_or_android", "macos", "windows"):
    from cpu_features.x86 cimport GetX86Info, X86Info
    cdef X86Info cpu_info = GetX86Info()
    _SSE2_BUILD_SUPPORT   = SSE2_BUILD_SUPPORT
    _SSE2_RUNTIME_SUPPORT = SSE2_BUILD_SUPPORT and cpu_info.features.sse2 != 0
    _SSE4_BUILD_SUPPORT   = SSE4_BUILD_SUPPORT
    _SSE4_RUNTIME_SUPPORT = SSE4_BUILD_SUPPORT and cpu_info.features.sse4_1 != 0
    _AVX2_BUILD_SUPPORT   = AVX2_BUILD_SUPPORT
    _AVX2_RUNTIME_SUPPORT = AVX2_BUILD_SUPPORT and cpu_info.features.avx2 != 0
ELIF TARGET_CPU == "arm":
    from cpu_features.arm cimport GetArmInfo, ArmInfo
    cdef ArmInfo arm_info = GetArmInfo()
    _NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
    _NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT and arm_info.features.neon != 0
ELIF TARGET_CPU == "aarch64":
    _NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
    _NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT # always runtime support on Aarch64

if _SSE2_RUNTIME_SUPPORT:
    _BEST_BACKEND = simd_backend.SSE2
elif _NEON_RUNTIME_SUPPORT:
    _BEST_BACKEND = simd_backend.NEON
else:
    _BEST_BACKEND = simd_backend.GENERIC

cdef enum simd_backend:
    NONE
    GENERIC
    SSE2
    SSE4
    AVX2
    NEON

IF SSE2_BUILD_SUPPORT:
    from pyopal._opal_sse2 cimport opalSearchDatabaseSSE2
IF SSE4_BUILD_SUPPORT:
    from pyopal._opal_sse4 cimport opalSearchDatabaseSSE4
IF AVX2_BUILD_SUPPORT:
    from pyopal._opal_avx2 cimport opalSearchDatabaseAVX2


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
) nogil


cdef class SearchResults:
    cdef OpalSearchResult** _results
    cdef size_t             _size

    def __cinit__(self):
        self._results = NULL
        self._size = 0

    def __dealloc__(self):
        PyMem_Free(self._results)

    def __len__(self):
        return self._size

    def __getitem__(self, ssize_t index):
        cdef ssize_t index_ = index
        if index_ < 0:
            index_ += self._size
        if index_ < 0 or index_ >= self._size:
            raise IndexError(index)
        
        cdef SearchResult result = SearchResult.__new__(SearchResult)
        result._result = self._results[index_]
        result._owner  = self
        return result


cdef class SearchResult:
    cdef SearchResults     _owner
    cdef OpalSearchResult* _result

    @property
    def score(self):
        assert self._result is not NULL
        return self._result.score if self._result.scoreSet else None

    @property
    def query_start(self):
        assert self._result is not NULL
        return None if self._result.startLocationQuery == -1 else self._result.startLocationQuery

    @property
    def query_end(self):
        assert self._result is not NULL
        return None if self._result.endLocationQuery == -1 else self._result.endLocationQuery

    @property
    def target_start(self):
        assert self._result is not NULL
        return None if self._result.startLocationTarget == -1 else self._result.startLocationTarget
        
    @property
    def target_end(self):
        assert self._result is not NULL
        return None if self._result.endLocationTarget == -1 else self._result.endLocationTarget

    @property
    def alignment(self):
        assert self._result is not NULL
        if self._result.alignment is NULL:
            return None
        return PyBytes_FromStringAndSize(<char*> self._result.alignment, self._result.alignmentLength)


cdef class Database:
    """A database of target sequences.
    """

    # cdef          size_t          _size
    # cdef          unsigned char** _targets
    # cdef          int*            _lengths
    cdef readonly tuple           sequences

    def __cinit__(self):
        self.sequences = tuple()
        # self._size      = 0
        # self._targets   = NULL
        # self._lengths   = NULL

    def __init__(self, *sequences):
        cdef size_t i
        cdef bytes  sequence
        
        self.sequences = tuple(sequences)
        # self._size     = len(self.sequences)
        # self._lengths  = <int*> PyMem_Malloc(self._size * sizeof(int))
        # self._targets  = <unsigned char**> PyMem_Malloc(self._size * sizeof(unsigned char*))
        #
        # if self._lengths is NULL or self._targets is NULL:
        #     raise MemoryError("Allocation failed")
        #
        # for i, sequence in enumerate(self.sequences):
        #     self._lengths[i] = len(sequence)
        #     self._targets[i] = <unsigned char*> PyBytes_AsString(sequence)

    def __dealloc__(self):
        pass
        # PyMem_Free(self._targets)
        # PyMem_Free(self._lengths)

    def search(self, bytes sequence, *, int gap_open = 3, int gap_extend = 1):
        cdef size_t          i
        cdef bytes           seq
        cdef int             retcode
        cdef searchfn_t      searchfn = NULL
        cdef SearchResults   results  = SearchResults.__new__(SearchResults)
        cdef ScoreMatrix     matrix   = ScoreMatrix.getBlosum50()

        # create encoding index
        cdef char            idx[255]
        memset(idx, 0, 255*sizeof(char))
        for i, letter in enumerate(matrix.getAlphabet()):
            idx[toupper(letter)] = idx[tolower(letter)] = i

        # prepare database
        cdef size_t          size    = len(self.sequences)
        cdef unsigned char** targets = <unsigned char**> PyMem_Malloc(len(self.sequences) * sizeof(unsigned char*))
        cdef int*            lengths = <int*> PyMem_Malloc(len(self.sequences) * sizeof(int))
        for i, seq in enumerate(self.sequences):
            lengths[i] = len(seq)
            targets[i] = <unsigned char*> PyMem_Malloc(len(seq)*sizeof(char))
            for j in range(len(seq)):
                targets[i][j] = idx[seq[j]]

        # Prepare query
        cdef size_t         length   = len(sequence)
        cdef unsigned char* query    = <unsigned char*> PyMem_Malloc(length*sizeof(char))
        for j in range(len(sequence)):
            query[j] = idx[sequence[j]]

        # Prepare array of results
        results._size    = size
        results._results = <OpalSearchResult**> PyMem_Malloc(size * sizeof(OpalSearchResult*))
        for i in range(size):
            results._results[i] = <OpalSearchResult*> PyMem_Malloc(sizeof(OpalSearchResult))
            opal.opalInitSearchResult(results._results[i])

        #if _AVX2_RUNTIME_SUPPORT:
        #    searchfn = opalSearchDatabaseAVX2
        #elif _SSE4_RUNTIME_SUPPORT:
        #    searchfn = opalSearchDatabaseSSE4
        #elif _SSE2_RUNTIME_SUPPORT:
        #else:
        #    raise RuntimeError("No supported SIMD backend available")

        searchfn = opalSearchDatabaseSSE2

        if searchfn is NULL:
            raise RuntimeError("Failed to detect SIMD backend")
    
        with nogil:
            retcode = searchfn(
                query,
                length,
                targets,
                size,
                lengths,
                gap_open,
                gap_extend,
                matrix.getMatrix(),
                matrix.getAlphabetLength(),
                results._results,
                opal.OPAL_SEARCH_ALIGNMENT,
                opal.OPAL_MODE_SW,
                opal.OPAL_OVERFLOW_BUCKETS,
            )

        if retcode != 0:
            raise RuntimeError("oh no!")

        return results