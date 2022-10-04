

cdef extern from "opal.h" nogil:
    
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