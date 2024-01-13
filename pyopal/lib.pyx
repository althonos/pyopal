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

from libc.string cimport memset, memcpy
from libc.stdint cimport uint32_t
from libc.limits cimport UCHAR_MAX
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from shared_mutex cimport shared_mutex

from cpython cimport Py_INCREF
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ, PyBUF_WRITE
from cpython.list cimport PyList_New, PyList_SET_ITEM
from cpython.bytes cimport PyBytes_AsString, PyBytes_FromStringAndSize
from cpython.mem cimport PyMem_Calloc, PyMem_Realloc, PyMem_Free
from cpython.memoryview cimport PyMemoryView_FromMemory
from cpython.unicode cimport (
    PyUnicode_KIND,
    PyUnicode_DATA,
    PyUnicode_READ,
    PyUnicode_GET_LENGTH,
    PyUnicode_1BYTE_KIND
)

cimport opal
cimport opal.score_matrix
from opal cimport OpalSearchResult

if NEON_BUILD_SUPPORT:
    from pyopal.platform.neon cimport opalSearchDatabaseNEON
if SSE2_BUILD_SUPPORT:
    from pyopal.platform.sse2 cimport opalSearchDatabaseSSE2
if SSE4_BUILD_SUPPORT:
    from pyopal.platform.sse4 cimport opalSearchDatabaseSSE4
if AVX2_BUILD_SUPPORT:
    from pyopal.platform.avx2 cimport opalSearchDatabaseAVX2

cdef extern from "<cctype>" namespace "std" nogil:
    cdef int toupper(int ch)
    cdef int tolower(int ch)
    cdef bint isalpha(int cht)

cdef extern from "<algorithm>" namespace "std" nogil:
    cdef void reverse[T](T, T)
    cdef int  count[T, V](T, T, V)

cdef extern from * nogil:
    """
    template<typename T>
    shared_ptr<T> pyshared(T* obj) {
        return shared_ptr<T>(obj, PyMem_Free);
    }
    """
    shared_ptr[T] pyshared[T](T* obj)

# --- Python imports -----------------------------------------------------------

import operator
from string import ascii_lowercase

include "_version.py"
include "_patch.pxi"

# --- Constants ----------------------------------------------------------------

cdef dict _OPAL_SEARCH_MODES = {
    "score": opal.OPAL_SEARCH_SCORE,
    "end": opal.OPAL_SEARCH_SCORE_END,
    "full": opal.OPAL_SEARCH_ALIGNMENT,
}

cdef dict _OPAL_SEARCH_RESULTS = {
    "score": ScoreResult,
    "end": EndResult,
    "full": FullResult,
}

cdef dict _OPAL_OVERFLOW_MODES = {
    "simple": opal.OPAL_OVERFLOW_SIMPLE,
    "buckets": opal.OPAL_OVERFLOW_BUCKETS,
}

cdef dict _OPAL_ALGORITHMS = {
    "nw": opal.OPAL_MODE_NW,
    "hw": opal.OPAL_MODE_HW,
    "ov": opal.OPAL_MODE_OV,
    "sw": opal.OPAL_MODE_SW,
}

cdef dict _OPAL_ALIGNMENT_OPERATION = {
    'M': opal.OPAL_ALIGN_MATCH,
    'D': opal.OPAL_ALIGN_DEL,
    'I': opal.OPAL_ALIGN_INS,
    'X': opal.OPAL_ALIGN_MISMATCH,
}

# --- Runtime CPU detection ----------------------------------------------------

if TARGET_SYSTEM == "windows":
    import cpuinfo
    _HOST_CPU = cpuinfo.get_cpu_info()
    _HOST_FEATURES = _HOST_CPU["flags"]
else:
    import archspec.cpu
    _HOST_CPU             = archspec.cpu.host()
    _HOST_FEATURES        = _HOST_CPU.features

_SSE2_BUILD_SUPPORT   = SSE2_BUILD_SUPPORT
_SSE4_BUILD_SUPPORT   = SSE4_BUILD_SUPPORT
_AVX2_BUILD_SUPPORT   = AVX2_BUILD_SUPPORT
_NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
_SSE2_RUNTIME_SUPPORT = SSE2_BUILD_SUPPORT and "sse2" in _HOST_FEATURES
_SSE4_RUNTIME_SUPPORT = SSE4_BUILD_SUPPORT and "sse4_1" in _HOST_FEATURES
_AVX2_RUNTIME_SUPPORT = AVX2_BUILD_SUPPORT and "avx2" in _HOST_FEATURES
_NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT and "neon" in _HOST_FEATURES

# NOTE(@althonos): NEON is always supported on Aarch64 so we should only check
#                  that the extension was built with NEON support.
if TARGET_CPU == "aarch64":
    _NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT

# NOTE(@althonos): SSE2 is always supported on x86-64 so we should only check
#                  that the extension was built with SSE2 support.
if TARGET_CPU == "x86_64":
    _SSE2_RUNTIME_SUPPORT = SSE2_BUILD_SUPPORT


# --- Exclusive read/write lock ------------------------------------------------

cdef class SharedMutex:

    def __cinit__(self):
        self.read = ReadLock(self)
        self.write = WriteLock(self)


cdef class ReadLock:

    def __cinit__(self, SharedMutex owner):
        self.owner = owner

    def __enter__(self):
        self.owner.mutex.lock_shared()

    def __exit__(self, exc_type, exc_value, traceback):
        self.owner.mutex.unlock_shared()


cdef class WriteLock:

    def __cinit__(self, SharedMutex owner):
        self.owner = owner

    def __enter__(self):
        self.owner.mutex.lock()

    def __exit__(self, exc_type, exc_value, traceback):
        self.owner.mutex.unlock()


# --- Python classes -----------------------------------------------------------

cdef class Alphabet:
    """A class for ordinal encoding of sequences.
    """

    def __init__(self, str letters = "ARNDCQEGHILKMFPSTWYVBZX*"):
        if len(letters) != len(set(letters)):
            raise ValueError("duplicate symbols in alphabet letters")
        if any(x != '*' and not x.isupper() for x in letters):
            raise ValueError("alphabet must only contain uppercase characters or wildcard")

        self.letters = letters
        self._unknown = letters.find('*')

        memset(&self._ahash[0], self._unknown, UCHAR_MAX * sizeof(char))
        for i, x in enumerate(letters.encode('ascii')):
            self._ahash[x] = i

    def __len__(self):
        return len(self.letters)

    def __contains__(self, object item):
        return item in self.letters

    def __getitem__(self, ssize_t index):
        return self.letters[index]

    def __reduce__(self):
        return type(self), (self.letters,)

    def __repr__(self):
        return f"{type(self).__name__}({self.letters!r})"

    def __str__(self):
        return self.letters

    def __eq__(self, object item):
        if isinstance(item, str):
            return self.letters == item
        elif isinstance(item, Alphabet):
            return self.letters == item.letters
        else:
            return False

    cdef void encode_raw(self, const unsigned char[:] sequence, digit_t[:] encoded):
        cdef ssize_t       i
        cdef char          code
        cdef unsigned char letter

        for i in range(sequence.shape[0]):
            letter = sequence[i]
            if not isalpha(letter):
                raise ValueError(f"character outside ASCII range: {letter!r}")
            code = self._ahash[<unsigned char> letter]
            if code < 0:
                raise ValueError(f"non-alphabet character in sequence: {chr(letter)!r}")
            encoded[i] = code

    cdef void decode_raw(self, const digit_t[:] encoded, unsigned char[:] sequence):
        cdef digit_t code
        cdef size_t  length = len(self.letters)
        for i in range(encoded.shape[0]):
            code = encoded[i]
            if code >= length:
                raise ValueError(f"invalid index in encoded sequence: {code!r}")
            sequence[i] = ord(self.letters[code])

    # ---

    cpdef bytes encode(self, object sequence):
        r"""Encode a sequence to an ordinal-encoded sequence using the alphabet.

        Arguments:
            sequence (`str` or byte-like object): The sequence to encode.

        Raises:
            `ValueError`: When the sequence contains invalid characters, or
            unknown sequence characters while the alphabet contains no
            wildcard character.

        Example:
            >>> alphabet = Alphabet("ACGT")
            >>> alphabet.encode("GATACA")
            b'\x02\x00\x03\x00\x01\x00'

        """
        cdef bytearray encoded = bytearray(len(sequence))
        if isinstance(sequence, str):
            sequence = sequence.encode('ascii')
        self.encode_raw(sequence, encoded)
        return bytes(encoded)

    cpdef str decode(self, object encoded):
        r"""Decode an ordinal-encoded sequence using the alphabet.

        Arguments:
            sequence (byte-like object): The sequence to decode.

        Raises:
            `ValueError`: When the sequence contains invalid indices.

        Example:
            >>> alphabet = Alphabet("ACGT")
            >>> alphabet.decode(bytearray([2, 0, 3, 0, 1, 0]))
            'GATACA'

        """
        cdef bytearray decoded = bytearray(len(encoded))
        self.decode_raw(encoded, decoded)
        return decoded.decode('ascii')


cdef class ScoreMatrix:
    """A score matrix for comparing character matches in sequences.

    Attributes:
        alphabet (`~pyopal.Alphabet`): The alphabet object storing the
            letters corresponding to column and row indices of the
            scoring matrix.

    """

    @classmethod
    def aa(cls):
        """Create a default amino-acid scoring matrix (BLOSUM50).
        """
        cdef int           i
        cdef char          unknown
        cdef unsigned char letter
        cdef size_t        length
        cdef ScoreMatrix   matrix  = ScoreMatrix.__new__(ScoreMatrix)

        matrix._mx = opal.score_matrix.ScoreMatrix.getBlosum50()

        length = matrix._mx.getAlphabetLength()
        matrix._shape[0] = matrix._shape[1] = length

        letters = PyBytes_FromStringAndSize(<char*> matrix._mx.getAlphabet(), length)
        matrix.alphabet = Alphabet(letters.decode('ascii'))
        return matrix

    def __init__(self, object alphabet not None, object matrix not None):
        """Create a new score matrix from the given alphabet and scores.

        Arguments:
            alphabet (`str` or `~pyopal.Alphabet`): The alphabet of the
                similarity matrix.
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
        cdef char          unknown
        cdef int           length    = len(matrix)
        cdef vector[uchar] abcvector = vector[uchar](length, 0)
        cdef vector[int]   scores    = vector[int](length*length, 0)

        if isinstance(alphabet, Alphabet):
            self.alphabet = alphabet
        elif isinstance(alphabet, str):
            self.alphabet = Alphabet(alphabet)
        else:
            raise TypeError(f"expected str or Alphabet, found {type(alphabet).__name__!r}")

        if len(self.alphabet) != length:
            raise ValueError("Alphabet must have the same length as matrix")
        if not all(len(row) == length for row in matrix):
            raise ValueError("`matrix` must be a square matrix")

        # FIXME: may be required implicitly by SIMD implementations
        # if length > 32:
        #     raise ValueError(f"Cannot use alphabet of more than 32 symbols: {alphabet!r}")

        # copy the alphabet for Opal
        for i, letter in enumerate(self.alphabet.letters):
            abcvector[i] = ord(letter)
        # copy the scores
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                scores[i*length+j] = value
        # create the Opal matrix
        self._mx = opal.score_matrix.ScoreMatrix(abcvector, scores)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.alphabet!r}, {self.matrix!r})"

    def __reduce__(self):
        return (type(self), (self.alphabet, self.matrix))

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        if flags & PyBUF_FORMAT:
            buffer.format = b"i"
        else:
            buffer.format = NULL
        buffer.buf = self._mx.getMatrix()
        buffer.internal = NULL
        buffer.itemsize = sizeof(int)
        buffer.len = self._shape[0] * self._shape[1] * sizeof(int)
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = <Py_ssize_t*> &self._shape
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def matrix(self):
        """`list` of `list` of `int`: The score matrix.
        """
        cdef int        i
        cdef int        j
        cdef int        length = self._mx.getAlphabetLength()
        cdef const int* scores = self._mx.getMatrix()

        return [
            [ scores[i*length + j] for j in range(length) ]
            for i in range(length)
        ]


cdef class ScoreResult:
    """The results of a search in ``score`` mode.
    """

    def __cinit__(self):
        self._target_index = -1
        self._result.scoreSet = 0
        self._result.startLocationQuery = -1
        self._result.startLocationTarget = -1
        self._result.endLocationQuery = -1
        self._result.endLocationTarget = -1
        self._result.alignmentLength = 0
        self._result.alignment = NULL

    def __dealloc__(self):
        PyMem_Free(self._result.alignment)

    def __init__(self, size_t target_index, int score):
        self._target_index = target_index
        self._result.score = score
        self._result.scoreSet = True

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.target_index}, score={self.score!r})"

    @property
    def target_index(self):
        """`int`: The index of the target in the database.
        """
        assert self._target_index >= 0
        return self._target_index

    @property
    def score(self):
        """`int`: The score of the alignment.
        """
        assert self._result.scoreSet
        return self._result.score


cdef class EndResult(ScoreResult):
    """The results of a search in ``end`` mode.
    """

    def __init__(
        self,
        size_t target_index,
        int score,
        int query_end,
        int target_end,
    ):
        super().__init__(target_index, score)
        self._result.endLocationQuery = query_end
        self._result.endLocationTarget = target_end

    def __repr__(self):
        cdef str ty = type(self).__name__
        return (
            f"{ty}({self.target_index}, "
            f"score={self.score!r}, "
            f"query_end={self.query_end!r}, "
            f"target_end={self.target_end!r})"
        )

    @property
    def query_end(self):
        """`int`: The coordinate where the alignment ends in the query.
        """
        assert self._result.endLocationQuery >= 0
        return self._result.endLocationQuery

    @property
    def target_end(self):
        """`int`: The coordinate where the alignment ends in the target.
        """
        assert self._result.endLocationTarget >= 0
        return self._result.endLocationTarget



cdef class FullResult(EndResult):
    """The results of a search in ``full`` mode.
    """

    def __cinit__(self):
        self._query_length = -1
        self._target_length = -1

    def __init__(
        self,
        size_t target_index,
        int score,
        int query_end,
        int target_end,
        int query_start,
        int target_start,
        int query_length,
        int target_length,
        str alignment not None,
    ):
        super().__init__(target_index, score, query_end, target_end)
        self._query_length = query_length
        self._target_length = target_length
        self._result.startLocationQuery = query_start
        self._result.startLocationTarget = target_start
        self._result.alignmentLength = len(alignment)
        self._result.alignment = <unsigned char*> PyMem_Realloc(self._result.alignment, self._result.alignmentLength * sizeof(unsigned char))
        for i, x in enumerate(alignment):
            self._result.alignment[i] = _OPAL_ALIGNMENT_OPERATION[x]

    def __repr__(self):
        cdef str ty = type(self).__name__
        return (
            f"{ty}({self.target_index}, "
            f"score={self.score!r}, "
            f"query_end={self.query_end!r}, "
            f"target_end={self.target_end!r}, "
            f"query_start={self.query_start!r}, "
            f"target_start={self.target_start!r}, "
            f"target_start={self.target_start!r}, "
            f"query_length={self.query_length!r}, "
            f"target_length={self.target_length!r}, "
            f"alignment={self.alignment!r})"
        )

    @property
    def query_start(self):
        """`int`: The coordinate where the alignment starts in the query.
        """
        assert self._result.startLocationQuery >= 0
        return self._result.startLocationQuery

    @property
    def target_start(self):
        """`int`: The coordinate where the alignment starts in the target.
        """
        assert self._result.startLocationTarget >= 0
        return self._result.startLocationTarget

    @property
    def query_length(self):
        """`int`: The complete length of the query sequence.

        .. versionadded:: 0.2.0

        """
        assert self._query_length >= 0
        return self._query_length

    @property
    def target_length(self):
        """`int`: The complete length of the target sequence.

        .. versionadded:: 0.2.0

        """
        assert self._target_length >= 0
        return self._target_length

    @property
    def alignment(self):
        """`str`: A string used by Opal to encode alignments.
        """
        assert self._result.alignmentLength >= 0

        cdef bytearray        ali
        cdef unsigned char[:] view
        cdef Py_UCS4[4]       symbols = ['M', 'D', 'I', 'X']

        ali = bytearray(self._result.alignmentLength)
        view = ali
        for i in range(self._result.alignmentLength):
            view[i] = symbols[self._result.alignment[i]]
        return ali.decode('ascii')

    cpdef str cigar(self):
        """Create a CIGAR string representing the alignment.

        Returns:
            `str`: A CIGAR string in SAM format describing the alignment.

        Example:
            >>> db = Database(["AACCGCTG"])
            >>> hit = db.search("ACCTCG", mode="full", algorithm="nw")[0]
            >>> hit.cigar()
            '1D5M1D1M'

        """
        cdef ssize_t       i
        cdef unsigned char symbol
        cdef unsigned char current
        cdef size_t        count
        cdef Py_UCS4[3]    symbols = ['M', 'I', 'D']
        cdef list          chunks  = []

        if self._result.alignmentLength == 0 or self._result.alignment is None:
            return None

        count = 0
        current = self._result.alignment[0]
        for i in range(self._result.alignmentLength):
            symbol = self._result.alignment[i] % 3
            if symbol == current:
                count += 1
            else:
                chunks.append(str(count))
                chunks.append(symbols[current])
                current = symbol
                count = 1
        chunks.append(str(count))
        chunks.append(symbols[current])

        return "".join(chunks)

    cpdef float identity(self):
        """Compute the identity of the alignment.

        Returns:
            `float`: The identity of the alignment as a fraction
            (between *0* and *1*).

        """
        assert self._result.alignment is not NULL
        cdef size_t         length = self._result.alignmentLength
        cdef unsigned char* ali    = self._result.alignment
        cdef int matches    = count(&ali[0], &ali[length], opal.OPAL_ALIGN_MATCH)
        cdef int mismatches = count(&ali[0], &ali[length], opal.OPAL_ALIGN_MISMATCH)
        return (<float> matches) / (<float> (matches + mismatches))

    cpdef float coverage(self, str reference="query"):
        """Compute the coverage of the alignment.

        Arguments:
            reference (`str`): The reference sequence to take to compute
                the coverage: either ``query`` or ``target``.

        Returns:
            `float`: The coverage of the alignment against the
            reference, as a fraction (between *0* and *1*).

        Example:
            If we align the test sequences from the Opal dataset with
            Needleman-Wunsch, we get the following alignment::

                T: AACCGCTG (0 - 7)
                Q: -ACCTC-G (0 - 5)

            The query coverage will be 100%, but the target coverage will
            only be 87.5%, since the first letter of the target is not
            covered by the alignment::

                >>> db = Database(["AACCGCTG"])
                >>> hit = db.search("ACCTCG", mode="full", algorithm="nw")[0]
                >>> hit.coverage("query")
                1.0
                >>> hit.coverage("target")
                0.875

        .. versionadded:: 0.2.0

        """
        assert self._result.alignment is not NULL

        cdef ssize_t i
        cdef ssize_t length
        cdef ssize_t reflength
        cdef char    operation

        # compute alignment and total lengths on reference sequence
        if reference == "query":
            reflength = self._query_length
            length = self._result.endLocationQuery + 1 - self._result.startLocationQuery
            operation = opal.OPAL_ALIGN_DEL
        elif reference == "target":
            reflength = self._target_length
            length = self._result.endLocationTarget + 1 - self._result.startLocationTarget
            operation = opal.OPAL_ALIGN_INS
        else:
            raise ValueError(f"Invalid coverage reference: {reference!r}")

        # trim alignment sides if they correspond to a gap in the reference
        for i in range(self._result.alignmentLength):
            if self._result.alignment[i] == operation:
                length -= 1
            else:
                break
        for i in reversed(range(self._result.alignmentLength)):
            if self._result.alignment[i] == operation:
                length -= 1
            else:
                break

        # compute the final coverage
        return 0.0 if length < 0 else (<float> length) / reflength


cdef class Database:
    """A database of target sequences.

    Like many biological sequence analysis tools, Opal encodes sequences
    with an alphabet for faster indexing of matrices. Sequences inserted in
    a database are stored in encoded format using the alphabet given on
    instantiation.

    Attributes:
        alphabet (`~pyopal.Alphabet`): The alphabet object used for encoding
            the sequences stored in the sequence database.

    """

    _DEFAULT_SCORE_MATRIX = ScoreMatrix.aa()

    # --- Magic methods --------------------------------------------------------

    def __cinit__(self):
        self.lock       = SharedMutex()
        self._sequences = vector[seq_t]()
        self._lengths   = vector[int]()
        self._search    = NULL

    def __init__(self, object sequences=(), Alphabet alphabet = None):
        # record the alphabet
        self.alphabet = self._DEFAULT_SCORE_MATRIX.alphabet if alphabet is None else alphabet

        # reset the collection if `__init__` is called more than once
        self.clear()
        # add the sequences to the database
        self.extend(sequences)

        # Select the best available SIMD backend
        if SSE2_BUILD_SUPPORT and _SSE2_RUNTIME_SUPPORT:
            self._search = opalSearchDatabaseSSE2
        if SSE4_BUILD_SUPPORT and _SSE4_RUNTIME_SUPPORT:
            self._search = opalSearchDatabaseSSE4
        if AVX2_BUILD_SUPPORT and _AVX2_RUNTIME_SUPPORT:
            self._search = opalSearchDatabaseAVX2
        if NEON_BUILD_SUPPORT and _NEON_RUNTIME_SUPPORT:
            self._search = opalSearchDatabaseNEON
        if self._search is NULL:
            raise RuntimeError("No supported SIMD backend available")

    def __reduce__(self):
        return (type(self), ((), self.alphabet), None, iter(self))

    # --- Encoding -------------------------------------------------------------

    cdef seq_t _encode(self, object sequence) except *:  
        cdef size_t   length  = len(sequence)
        cdef bytes    encoded = self.alphabet.encode(sequence)
        cdef char*    indices = encoded
        cdef digit_t* dst     = NULL

        dst = <digit_t*> PyMem_Calloc(length, sizeof(digit_t))
        if dst is NULL:
            raise MemoryError("Failed to allocate sequence data")

        memcpy(<void*> &dst[0], <void*> indices, length * sizeof(digit_t))
        return pyshared(dst)

    cdef str _decode(self, seq_t encoded, int length) except *:
        view = PyMemoryView_FromMemory(<char*> encoded.get(), length, PyBUF_READ)
        return self.alphabet.decode(view)

    # --- Sequence interface ---------------------------------------------------

    def __len__(self):
        with self.lock.read:
            return self._sequences.size()

    def __getitem__(self, ssize_t index):
        cdef ssize_t i
        cdef size_t  size
        cdef object  view
        cdef ssize_t index_   = index

        with self.lock.read:
            size = self._sequences.size()

            if index_ < 0:
                index_ += size
            if index_ < 0 or (<size_t> index_) >= size:
                raise IndexError(index)

            return self._decode(self._sequences[index_], self._lengths[index_])

    def __setitem__(self, ssize_t index, object sequence):
        cdef size_t  size
        cdef ssize_t index_  = index
        cdef int     length  = len(sequence)
        cdef seq_t   encoded = self._encode(sequence)

        with self.lock.write:
            size = self._sequences.size()

            if index_ < 0:
                index_ += size
            if index_ < 0 or (<size_t> index_) >= size:
                raise IndexError(index)

            self._sequences[index_] = encoded
            self._lengths[index_] = length
            self._pointers[index_] = encoded.get()

    def __delitem__(self, ssize_t index):
        cdef int     length
        cdef seq_t   encoded
        cdef size_t  size
        cdef ssize_t index_  = index

        with self.lock.write:
            size = self._sequences.size()

            if index_ < 0:
                index_ += size
            if index_ < 0 or (<size_t> index_) >= size:
                raise IndexError(index)

            self._sequences.erase(self._sequences.begin() + index_)
            self._lengths.erase(self._lengths.begin() + index_)
            self._pointers.erase(self._pointers.begin() + index_)

    cpdef void clear(self) except *:
        """Remove all sequences from the database.
        """
        cdef size_t i

        with self.lock.write:
            self._sequences.clear()
            self._lengths.clear()
            self._pointers.clear()

    cpdef void extend(self, object sequences) except *:
        """Extend the database by adding sequences from an iterable.

        Example:
            >>> db = pyopal.Database(["ATGC"])
            >>> db.extend(["TTCA", "AAAA", "GGTG"])
            >>> list(db)
            ['ATGC', 'TTCA', 'AAAA', 'GGTG']

        """
        cdef size_t size
        cdef size_t hint = operator.length_hint(sequences)

        # attempt to reserve space in advance for the new sequences
        with self.lock.write:
            size = self._sequences.size()
            if hint > 0:
                self._sequences.reserve(hint + size)
                self._lengths.reserve(hint + size)
                self._pointers.reserve(hint + size)

        # append sequences in order
        for sequence in sequences:
            self.append(sequence)

    cpdef void append(self, object sequence) except *:
        """Append a single sequence at the end of the database.

        Arguments:
            sequence (`str` or byte-like object): The new sequence.

        Hint:
            When inserting several sequences in the database, consider
            using the `Database.extend` method instead so that the
            internal buffers can reserve space just once for every
            new sequence.

        Example:
            >>> db = pyopal.Database(["ATGC", "TTCA"])
            >>> db.append("AAAA")
            >>> list(db)
            ['ATGC', 'TTCA', 'AAAA']

        """
        cdef int   length  = len(sequence)
        cdef seq_t encoded = self._encode(sequence)

        with self.lock.write:
            self._sequences.push_back(encoded)
            self._lengths.push_back(length)
            self._pointers.push_back(encoded.get())

    cpdef void reverse(self) except *:
        """Reverse the database, in place.

        Example:
            >>> db = pyopal.Database(['ATGC', 'TTGC', 'CTGC'])
            >>> db.reverse()
            >>> list(db)
            ['CTGC', 'TTGC', 'ATGC']

        """
        with self.lock.write:
            reverse(self._sequences.begin(), self._sequences.end())
            reverse(self._lengths.begin(), self._lengths.end())
            reverse(self._pointers.begin(), self._pointers.end())

    cpdef void insert(self, ssize_t index, object sequence):
        """Insert a sequence in the database at a given position.

        Arguments:
            index (`int`): The index where to insert the new sequence.
            sequence (`str` or byte-like object): The new sequence.

        Note:
            If the insertion index is out of bounds, the insertion will
            happen at either end of the database::

                >>> db = pyopal.Database(["ATGC", "TTGC", "CTGC"])
                >>> db.insert(-100, "TTTT")
                >>> db.insert(100, "AAAA")
                >>> list(db)
                ['TTTT', 'ATGC', 'TTGC', 'CTGC', 'AAAA']

        """
        cdef size_t  size
        cdef int     length  = len(sequence)
        cdef ssize_t index_  = index
        cdef seq_t   encoded = self._encode(sequence)

        with self.lock.write:
            size = self._sequences.size()

            if index_ < 0:
                index_ += size

            if index_ < 0:
                index_ = 0
            elif (<size_t> index_) >= size:
                index_ = size

            self._sequences.insert(self._sequences.begin() + index_, encoded)
            self._lengths.insert(self._lengths.begin() + index_, length)
            self._pointers.insert(self._pointers.begin() + index_, encoded.get())


    # --- Subset ---------------------------------------------------------------

    cpdef Database mask(self, object bitmask):
        """Extract a subset of the database where the bitmask is `True`.

        Arguments:
            bitmask (`collections.abc.Sequence` of `bool`): A sequence of
                `bool` objects with the same length as the database.

        Raises:
            `IndexError`: When the bitmask has a different dimension.

        Example:
            >>> db = pyopal.Database(['AAAA', 'CCCC', 'KKKK', 'FFFF'])
            >>> list(db.mask([True, False, False, True]))
            ['AAAA', 'FFFF']

        .. versionadded:: 0.3.0

        """
        cdef ssize_t  i
        cdef int      length
        cdef seq_t    seq
        cdef Database subdb

        if len(bitmask) != len(self):
            raise IndexError(bitmask)

        subdb = Database.__new__(Database)
        subdb.alphabet = self.alphabet
        subdb._search = self._search

        for i, b in enumerate(bitmask):
            if b:
                subdb._sequences.push_back(self._sequences[i])
                subdb._lengths.push_back(self._lengths[i])
                subdb._pointers.push_back(subdb._sequences[i].get())

        return subdb

    cpdef Database extract(self, object indices):
        """Extract a subset of the database using the given indices.

        Arguments:
            indices (`collections.abc.Sequence` of `int`): A sequence of
                `int` objects to use to index the database.

        Raises:
            `IndexError`: When ``indices`` contains an invalid index.

        Example:
            >>> db = pyopal.Database(['AAAA', 'CCCC', 'KKKK', 'FFFF'])
            >>> list(db.extract([2, 0]))
            ['KKKK', 'AAAA']

        .. versionadded:: 0.3.0

        """
        cdef ssize_t  index
        cdef int      length
        cdef seq_t    seq
        cdef Database subdb

        subdb = Database.__new__(Database)
        subdb.alphabet = self.alphabet
        subdb._search = self._search
        subdb._sequences.reserve(len(indices))
        subdb._pointers.reserve(len(indices))
        subdb._lengths.reserve(len(indices))

        for index in indices:
            if index < 0 or index >= len(self):
                raise IndexError(index)
            subdb._sequences.push_back(self._sequences[index])
            subdb._lengths.push_back(self._lengths[index])
            subdb._pointers.push_back(subdb._sequences[index].get())

        return subdb

    # --- Opal search ----------------------------------------------------------

    def search(
        self,
        object sequence,
        ScoreMatrix score_matrix = None,
        *,
        int gap_open = 3,
        int gap_extend = 1,
        str mode = "score",
        str overflow = "buckets",
        str algorithm = "sw",
    ):
        """Search the database with a query sequence.

        Arguments:
            sequence (`str` or byte-like object): The sequence to query the
                database with.

        Keyword Arguments:
            gap_open(`int`): The gap opening penalty :math:`G`.
            gap_extend (`int`): The gap extension penalty :math:`E`.
            mode (`str`): The search mode to use for querying the database:
                ``score`` to only report scores for each hit (default),
                ``end`` to report scores and end coordinates for each
                hit (slower), ``full`` to report scores, coordinates and
                alignment for each hit (slowest).
            overflow (`str`): The strategy to use when a sequence score
                overflows in the comparison pipeline: ``simple`` computes
                scores with 8-bit range first then recomputes with 16-bit
                range (and then 32-bit) the sequences that overflowed;
                ``buckets`` to divide the targets in buckets, and switch
                to larger score ranges within a bucket when the first
                overflow is detected.
            algorithm (`str`): The alignment algorithm to use: ``nw``
                for global Needleman-Wunsch alignment, ``hw`` for semi-global
                alignment without penalization of gaps on query edges, ``ov``
                for semi-global alignment without penalization of gaps on
                query or target edges, and ``sw`` for local Smith-Waterman
                alignment.

        Hit:
            A gap of length :math:`N` will receive a penalty of
            :math:`E + (N - 1)G`.

        Returns:
            `list` of `pyopal.ScoreResult`: A list containing one
                `ScoreResult` object for each target sequence in the
                database. The actual type depends on the requested
                ``mode``: it will be a `ScoreResult` for mode ``score``,
                `EndResult` for mode ``end`` and `FullResult` for mode
                ``full``.

        Raises:
            `ValueError`: When ``sequence`` contains invalid characters
                with respect to the alphabet of the database scoring
                matrix.
            `OverflowError`: When the score computed by Opal for a
                sequence overflows or underflows the limit values for
                the SIMD backend.
            `RuntimeError`: When attempting to call `Database.search` on
                a platform with no supported SIMD backend.

        """
        assert self.alphabet is not None
        assert self._search is not NULL

        cdef int                       _mode
        cdef int                       _overflow
        cdef int                       _algo
        cdef size_t                    i
        cdef int                       retcode
        cdef ScoreResult               result
        cdef FullResult                full_result
        cdef type                      result_type
        cdef list                      results
        cdef vector[OpalSearchResult*] results_raw
        cdef size_t                    size
        cdef seq_t                     query
        cdef int                       length      = len(sequence)

        # validate parameters
        if mode in _OPAL_SEARCH_MODES:
            _mode = _OPAL_SEARCH_MODES[mode]
            result_type = _OPAL_SEARCH_RESULTS[mode]
        else:
            raise ValueError(f"Invalid search mode: {mode!r}")
        if overflow in _OPAL_OVERFLOW_MODES:
            _overflow = _OPAL_OVERFLOW_MODES[overflow]
        else:
            raise ValueError(f"Invalid overflow mode: {overflow!r}")
        if algorithm in _OPAL_ALGORITHMS:
            _algo = _OPAL_ALGORITHMS[algorithm]
        else:
            raise ValueError(f"Invalid algorithm: {algorithm!r}")

        # check score matrix
        if score_matrix is None:
            score_matrix = self._DEFAULT_SCORE_MATRIX
        if score_matrix.alphabet != self.alphabet:
            raise ValueError("database and score matrix have different alphabets")

        # encode query
        query = self._encode(sequence)

        # search database
        with self.lock.read:
            size = self._sequences.size()
            # Prepare list of results
            res_array = PyMem_Calloc(size, sizeof(OpalSearchResult*))
            results_raw.reserve(size)
            results = PyList_New(size)
            for i in range(size):
                result = result_type.__new__(result_type)
                result._target_index = i
                Py_INCREF(result)
                PyList_SET_ITEM(results, i, result)
                results_raw.push_back(&result._result)
            # Run search pipeline in nogil mode
            with nogil:
                retcode = self._search(
                    query.get(),
                    length,
                    self._pointers.data(),
                    self._pointers.size(),
                    self._lengths.data(),
                    gap_open,
                    gap_extend,
                    score_matrix._mx.getMatrix(),
                    score_matrix._mx.getAlphabetLength(),
                    results_raw.data(),
                    _mode,
                    _algo,
                    _overflow,
                )

        # record query and target lengths if in full mode so that
        # the alignment coverage can be computed later
        if _mode == opal.OPAL_SEARCH_ALIGNMENT:
            for i in range(size):
                full_result = results[i]
                full_result._query_length = length
                full_result._target_length = self._lengths[i]

        # check the alignment worked and return results
        if retcode == opal.OPAL_ERR_NO_SIMD_SUPPORT:
            raise RuntimeError("No available SIMD on this platform")
        elif retcode == opal.OPAL_ERR_OVERFLOW:
            raise OverflowError("Overflow detected while computing alignment scores")
        elif retcode != 0:
            raise RuntimeError(f"Failed to run search Opal database (code={retcode})")
        return results







