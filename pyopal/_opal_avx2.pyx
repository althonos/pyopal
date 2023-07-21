# distutils: language = c++
# cython: language_level=3, linetrace=True, embedsignature=True, binding=True

cimport opal

include "_version.py"
include "_patch.pxi"

cdef int opalSearchDatabaseAVX2(
    unsigned char query[], 
    int queryLength, 
    unsigned char* db[], 
    int dbLength,
    int dbSeqLengths[], 
    int gapOpen, 
    int gapExt, 
    int* scoreMatrix,
    int alphabetLength, 
    opal.OpalSearchResult* results[],
    const int searchType, 
    int mode, 
    int overflowMethod
) nogil:
    return opal.opalSearchDatabase(
        query,
        queryLength, 
        db, 
        dbLength,
        dbSeqLengths, 
        gapOpen, 
        gapExt, 
        scoreMatrix,
        alphabetLength, 
        results,
        searchType, 
        mode, 
        overflowMethod
    )

