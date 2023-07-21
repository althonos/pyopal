

cdef extern from "opal.h" nogil:

    const int OPAL_ERR_OVERFLOW
    const int OPAL_ERR_NO_SIMD_SUPPORT
    const int OPAL_ERR_INVALID_MODE

    const int OPAL_MODE_NW
    const int OPAL_MODE_HW
    const int OPAL_MODE_OV
    const int OPAL_MODE_SW

    const int OPAL_OVERFLOW_SIMPLE
    const int OPAL_OVERFLOW_BUCKETS

    const int OPAL_SEARCH_SCORE
    const int OPAL_SEARCH_SCORE_END
    const int OPAL_SEARCH_ALIGNMENT

    const int OPAL_ALIGN_MATCH
    const int OPAL_ALIGN_DEL
    const int OPAL_ALIGN_INS
    const int OPAL_ALIGN_MISMATCH

    cdef struct OpalSearchResult:
        bint scoreSet
        int score
        int endLocationTarget
        int endLocationQuery
        int startLocationTarget
        int startLocationQuery
        unsigned char* alignment
        int alignmentLength

    cdef void opalInitSearchResult(OpalSearchResult* result)
    cdef bint opalSearchResultIsEmpty(const OpalSearchResult result)
    cdef void opalSearchResultSetScore(OpalSearchResult* result, int score)

    cdef int opalSearchDatabase(
        unsigned char query[],
        int queryLength,
        unsigned char* db[],
        int dbLength,
        int dbSeqLengths[],
        int gapOpen,
        int gapExt,
        int* scoreMatrix,
        int alphabetLength,
        OpalSearchResult* results[],
        const int searchType,
        int mode,
        int overflowMethod
    )

    cdef int opalSearchDatabaseCharSW(
        unsigned char query[],
        int queryLength,
        unsigned char** db,
        int dbLength,
        int dbSeqLengths[],
        int gapOpen,
        int gapExt,
        int* scoreMatrix,
        int alphabetLength,
        OpalSearchResult* results[]
    )