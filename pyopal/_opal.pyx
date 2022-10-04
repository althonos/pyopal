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
from libc.limits cimport UCHAR_MAX
from libcpp.vector cimport vector

cimport opal
cimport opal.score_matrix
from opal cimport OpalSearchResult

IF NEON_BUILD_SUPPORT:
    from pyopal._opal_neon cimport opalSearchDatabaseNEON
IF SSE2_BUILD_SUPPORT:
    from pyopal._opal_sse2 cimport opalSearchDatabaseSSE2
IF SSE4_BUILD_SUPPORT:
    from pyopal._opal_sse4 cimport opalSearchDatabaseSSE4
IF AVX2_BUILD_SUPPORT:
    from pyopal._opal_avx2 cimport opalSearchDatabaseAVX2

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


# --- Type definitions ---------------------------------------------------------

ctypedef unsigned char digit_t
ctypedef digit_t*      seq_t

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


# --- Sequence encoding --------------------------------------------------------

ctypedef fused bytelike:
    str
    object

cdef inline void encode(bytelike sequence, char* lut, seq_t* encoded, int* length) except *:

    cdef size_t                 i
    cdef char                   code
    cdef unsigned char          letter
    cdef const unsigned char[:] view

    if bytelike is str:
        raise NotImplementedError("Encoding from string not yet implemented")

    else:

        view       = sequence
        length[0]  = view.shape[0]
        encoded[0] = <seq_t> PyMem_Malloc(length[0] * sizeof(digit_t))
        if encoded[0] is NULL:
            raise MemoryError("Failed to allocate sequence data")

        for i in range(length[0]):
            letter = <unsigned char> view[i]
            code = lut[letter]
            if code < 0:
                raise ValueError(f"Invalid character in sequence: {chr(letter)}")
            encoded[0][i] = code

    print("encoded:", [encoded[0][i] for i in range(length[0])])


# --- Python classes -----------------------------------------------------------

cdef class ScoreMatrix:
    """A score matrix for comparing character matches in sequences.
    """
    cdef opal.score_matrix.ScoreMatrix _mx
    cdef char                          _ahash[UCHAR_MAX]

    @classmethod
    def aa(cls):
        """aa(cls)\n--

        Create a default amino-acid scoring matrix (BLOSUM50).

        """
        cdef size_t               i
        cdef unsigned char        letter
        cdef const unsigned char* alphabet

        cdef ScoreMatrix matrix = ScoreMatrix.__new__(ScoreMatrix)
        matrix._mx = opal.score_matrix.ScoreMatrix.getBlosum50()

        alphabet = matrix._mx.getAlphabet()

        for i in range(UCHAR_MAX):
            matrix._ahash[i] = -1
        for i in range(matrix._mx.getAlphabetLength()):
            letter = alphabet[i]
            matrix._ahash[toupper(letter)] = matrix._ahash[tolower(letter)] = i

        return matrix

    def __cinit__(self):
        memset(self._ahash, 0xFF, UCHAR_MAX*sizeof(char))

    def __init__(self, str alphabet not None, object matrix not None):
        """__init__(alphabet, matrix)\n--

        Create a new score matrix from the given alphabet and scores.

        Arguments:
            alphabet (`str`): The alphabet of the similarity matrix.
            matrix (`~numpy.typing.ArrayLike` of `int`): The scoring matrix,
                as a square matrix indexed by the alphabet characters.

        Example:
            Create a new similarity matrix using the HOXD70 scores by
            Chiaromonte, Yap and Miller (:pmid:`11928468`)::

                >>> matrix = ScoreMatrix(
                ...     "ATCG",
                ...     [[  91, -114,  -31, -123],
                ...      [-114,  100, -125,  -31],
                ...      [ -31, -125,  100, -114],
                ...      [-123,  -31, -114,   91]]
                ... )

            Create a new similarity matrix using one of the matrices from
            the `Bio.Align.substitution_matrices` module::

                >>> jones = Bio.Align.substitution_matrices.load('JONES')
                >>> matrix = ScoreMatrix(jones.alphabet, jones)

        """
        cdef int           i
        cdef int           j
        cdef object        row
        cdef int           value
        cdef str           letter
        cdef int           length    = len(matrix)
        cdef vector[uchar] abcvector = vector[uchar](length, 0)
        cdef vector[int]   scores    = vector[int](length*length, 0)

        if len(set(alphabet)) != len(alphabet):
            raise ValueError("Alphabet contains duplicate letters")
        if len(alphabet) != length:
            raise ValueError("Alphabet must have the same length as matrix")
        if not alphabet.isupper():
            raise ValueError("Alphabet must only contain uppercase letters")
        if not all(len(row) == length for row in matrix):
            raise ValueError("`matrix` must be a square matrix")

        # FIXME: may be required implicitly by SIMD implementations
        # if length > 32:
        #     raise ValueError(f"Cannot use alphabet of more than 32 symbols: {alphabet!r}")

        # copy the alphabet and create a lookup-table for encoding sequences
        for i, letter in enumerate(alphabet):
            j = ord(letter)
            abcvector[i] = j
            self._ahash[toupper(j)] = self._ahash[tolower(j)] = i

        # copy the scores
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                scores[i*length+j] = value

        # record the matrix
        self._mx = opal.score_matrix.ScoreMatrix(abcvector, scores)


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

    Like many biological sequence analysis tools, Opal encodes sequences
    with an alphabet for faster indexing of matrices. Sequences inserted in
    a  database are stored in encoded format using the alphabet of the
    `ScoreMatrix` given on instantiation.

    """

    cdef readonly ScoreMatrix   score_matrix
    cdef          vector[seq_t] _sequences
    cdef          vector[int]   _lengths

    # --- Magic methods --------------------------------------------------------

    def __cinit__(self):
        self._sequences = vector[seq_t]()
        self._lengths   = vector[int]()

    def __init__(self, object sequences=(), ScoreMatrix score_matrix=None):
        # reset the collection is `__init__` is called more than once
        self.clear()
        # record the score matrix
        self.score_matrix = score_matrix or ScoreMatrix.aa()
        # add the sequences to the database
        self.extend(sequences)

    def __dealloc__(self):
        cdef size_t i
        for i in range(self._sequences.size()):
            PyMem_Free(self._sequences[i])


    # --- Sequence interface ---------------------------------------------------

    def __len__(self):
        return self._sequences.size()

    cpdef void clear(self) except *:
        """clear(self)\n--

        Remove all sequences from the database.

        """
        cdef size_t i
        for i in range(self._sequences.size()):
            PyMem_Free(self._sequences[i])
        self._sequences.clear()
        self._lengths.clear()

    cpdef void extend(self, object sequences) except *:
        """extend(self, sequences)\n--

        Extend the database by adding sequences from an iterable.

        """
        for sequence in sequences:
            self.append(sequence)

    cpdef void append(self, bytelike sequence) except *:
        """append(self, sequence)\n--

        Append a single sequence at the end of the database.

        """
        cdef int   length
        cdef seq_t encoded

        encode(sequence, self.score_matrix._ahash, &encoded, &length)

        self._sequences.push_back(encoded)
        self._lengths.push_back(length)

    # --- Opal search ----------------------------------------------------------

    def search(self, bytelike sequence, *, int gap_open = 3, int gap_extend = 1):
        cdef size_t          i
        cdef int             retcode
        cdef int             length
        cdef size_t          size     = self._sequences.size()
        cdef seq_t           query    = NULL
        cdef searchfn_t      searchfn = NULL
        cdef SearchResults   results  = SearchResults.__new__(SearchResults)

        # Prepare query
        encode(sequence, self.score_matrix._ahash, &query, &length)

        # Prepare array of results
        results._size    = size
        results._results = <OpalSearchResult**> PyMem_Malloc(size * sizeof(OpalSearchResult*))
        for i in range(size):
            results._results[i] = <OpalSearchResult*> PyMem_Malloc(sizeof(OpalSearchResult))
            opal.opalInitSearchResult(results._results[i])

        # Select best available SIMD backend
        IF AVX2_BUILD_SUPPORT:
            if _AVX2_RUNTIME_SUPPORT and searchfn is NULL:
                searchfn = opalSearchDatabaseAVX2
        IF SSE4_BUILD_SUPPORT:
            if _SSE4_RUNTIME_SUPPORT and searchfn is NULL:
                searchfn = opalSearchDatabaseSSE4
        IF SSE2_BUILD_SUPPORT:
            if _SSE2_RUNTIME_SUPPORT and searchfn is NULL:
                searchfn = opalSearchDatabaseSSE2
        if searchfn is NULL:
            raise RuntimeError("No supported SIMD backend available")

        # Run search pipeline in nogil mode
        with nogil:
            retcode = searchfn(
                query,
                length,
                self._sequences.data(),
                self._sequences.size(),
                self._lengths.data(),
                gap_open,
                gap_extend,
                self.score_matrix._mx.getMatrix(),
                self.score_matrix._mx.getAlphabetLength(),
                results._results,
                opal.OPAL_SEARCH_ALIGNMENT,
                opal.OPAL_MODE_SW,
                opal.OPAL_OVERFLOW_BUCKETS,
            )

        # free allocated memory
        PyMem_Free(query)

        # return results or raise an exception on failure
        if retcode != 0:
            raise RuntimeError(f"Failed to run search Opal database (code={retcode})")
        return results







