

cdef extern from "opal.h" nogil:

    const int OPAL_ERR_OVERFLOW        = 1
    const int OPAL_ERR_NO_SIMD_SUPPORT = 2
    const int OPAL_ERR_INVALID_MODE    = 3

    const int OPAL_MODE_NW = 0
    const int OPAL_MODE_HW = 1
    const int OPAL_MODE_OV = 2
    const int OPAL_MODE_SW = 3

    const int OPAL_OVERFLOW_SIMPLE  = 0
    const int OPAL_OVERFLOW_BUCKETS = 1

    const int OPAL_SEARCH_SCORE     = 0
    const int OPAL_SEARCH_SCORE_END = 1
    const int OPAL_SEARCH_ALIGNMENT = 2

    const int OPAL_ALIGN_MATCH    = 0
    const int OPAL_ALIGN_DEL      = 1
    const int OPAL_ALIGN_INS      = 2
    const int OPAL_ALIGN_MISMATCH = 3
    
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