# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True

from cython.operator cimport dereference, preincrement

cimport libcpp.algorithm
from libc.stdint cimport int32_t, uint16_t, uint32_t, uint64_t
from libcpp cimport bool, nullptr
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr

cimport sword.chain
cimport sword.hash
cimport sword.kmers
cimport sword.reader
cimport sword.score_matrix
from sword.database_search cimport ChainEntry as _ChainEntry, ChainEntrySet as _ChainEntrySet, Indexes as _Indexes
from sword.kmers cimport Kmers as _Kmers
from sword.chain cimport ChainSet as _ChainSet, Chain as _Chain
from sword.reader cimport Reader as _Reader 
from sword.hash cimport Iterator as _HashIterator
from sword.score_matrix cimport ScoreMatrix as _ScoreMatrix, ScoreMatrixType as _ScoreMatrixType

import os


cdef class Kmers:
    cdef shared_ptr[_Kmers] _kmers

    def __init__(self, ScoreMatrix score_matrix, kmer_length = 3, score_threshold = 13):
        self._kmers = shared_ptr[_Kmers](
            sword.kmers.createKmers(kmer_length, score_threshold, score_matrix._score_matrix)
        )


cdef class ChainSet:
    cdef _ChainSet _chains

    def __len__(self):
        return self._chains.size()

    def __getitem__(self, ssize_t i):
        cdef ssize_t _i = i
        if _i < 0:
            _i += self._chains.size()
        if _i < 0 or _i >= self._chains.size():
            raise IndexError(i)
        return self._chains[i].get().name()


cdef class Reader:
    cdef shared_ptr[_Reader] _reader

    def __init__(self, path):
        cdef bytes _path = os.fsencode(path)
        self._reader = shared_ptr[_Reader](sword.reader.createReader(_path))

    cpdef ChainSet read(self):
        assert self._reader != nullptr
        cdef ChainSet s = ChainSet()
        self._reader.get().read_chains(s._chains, 0)
        return s


cdef class ScoreMatrix:
    cdef shared_ptr[_ScoreMatrix] _score_matrix

    def __init__(self, str name = "BLOSUM62", int32_t gap_open = 10, int32_t gap_extend = 1):
        cdef _ScoreMatrixType ty
        if name == "BLOSUM62":
            ty = sword.score_matrix.kBlosum62
        else:
            raise ValueError(f"unsupported score matrix: {name!r}")
        self._score_matrix = shared_ptr[_ScoreMatrix](
            sword.score_matrix.createScoreMatrix(
                ty,
                gap_open,
                gap_extend,
            )
        )

    @property
    def name(self):
        assert self._score_matrix != nullptr
        return self._score_matrix.get().scorerName().decode()


cdef extern from * nogil:
    """
    bool chainLengthKey(const std::unique_ptr<Chain>& left, const std::unique_ptr<Chain>& right) {
        return left->length() < right->length();
    }
    """
    bool chainLengthKey(const unique_ptr[_Chain]& left, const unique_ptr[_Chain]& right) noexcept


cpdef vector[uint32_t] preprocess_database(
    ChainSet database,
    size_t num_threads = 1,
    uint32_t max_short_length = 2000,
):
    # ported from: preprocDatabase in `database_search.cpp`
    cdef vector[uint32_t] dst
    cdef uint32_t         i

    # sort by length
    libcpp.algorithm.sort(database._chains.begin(), database._chains.end(), chainLengthKey)

    # split tasks between long and short
    cdef uint64_t short_total_length = 0
    cdef uint64_t long_total_length  = 0
    cdef uint32_t split              = 0

    for i in range(database._chains.size()):
        l = database._chains[i].get().length() 
        if l > max_short_length:
            if split == 0:
                split = i
            long_total_length += l
        else:
            short_total_length += l

    if short_total_length == 0:
        split = 0
    if long_total_length == 0:
        split = database._chains.size()

    # spread tasks across threads
    cdef uint64_t short_task_size = short_total_length / num_threads
    cdef uint64_t long_task_size = long_total_length / num_threads

    dst.reserve(2*num_threads + 1)
    dst.emplace_back(0)

    #
    cdef uint64_t total_length = 0

    for i in range(split):
        total_length += database._chains[i].get().length()
        if total_length > short_task_size:
            total_length = 0
            dst.emplace_back(i + 1)
    
    if dst.back() != split:
        dst.emplace_back(split)

    total_length = 0
    for i in range(split, database._chains.size()):
        total_length += database._chains[i].get().length()
        if total_length > long_task_size:
            total_length = 0
            dst.emplace_back(i + 1)

    if dst.back() != database._chains.size():
        dst.emplace_back(database._chains.size())

    return dst


cdef extern from * nogil:
    """
    bool chainEntryDataKey(const ChainEntry& left, const ChainEntry& right) {
        return left.data() > right.data();
    }
    """
    bool chainEntryDataKey(const _ChainEntry& left, const _ChainEntry& right) noexcept


cdef uint32_t    kProtBits   = 5
cdef uint32_t[6] kmerDelMask = [ 0, 0, 0, 0x7fff, 0xFFFFF, 0x1FFFFFF ]

cdef _ChainEntrySet score_chains(
    ChainSet queries,
    ChainSet database,
    size_t max_candidates,
    uint32_t database_start,
    uint32_t database_end,
    Kmers kmers,
):
    cdef uint32_t      kmer
    cdef _HashIterator begin
    cdef _HashIterator end
    
    cdef _ChainEntrySet entries_part = _ChainEntrySet(queries._chains.size())
    cdef vector[uint16_t] min_entry_score = vector[uint16_t](queries._chains.size())
    cdef vector[uint16_t] entries_found   = vector[uint16_t](queries._chains.size())

    cdef _ChainEntrySet dst = _ChainEntrySet(queries._chains.size())

    cdef uint32_t id_
    for i in range(queries._chains.size()):
        id_ = queries._chains[i].get().id()
        entries_found[i] = dst[id_].size()
        min_entry_score[i] = 65000 if entries_found[i] == 0 else dst[id_].back().data()

    cdef uint32_t kmer_length       = kmers._kmers.get().kmer_length()
    cdef uint32_t max_target_length = database._chains[database_end - 1].get().length()
    cdef size_t   groups            = 0

    cdef uint32_t kmer_offset       = kmer_length - 1
    cdef uint32_t del_mask          = kmerDelMask[kmer_length]

    cdef uint32_t max_scores_length = 100000 if kmer_length == 3 else 500000;

    cdef vector[uint16_t] scores = vector[uint16_t](max_scores_length)
    cdef vector[uint32_t] score_lengths = vector[uint32_t](queries._chains.size())
    cdef vector[uint32_t] score_starts = vector[uint32_t](queries._chains.size() + 1)
    cdef vector[uint16_t] max_score = vector[uint16_t](queries._chains.size())
    score_starts[0] = 0

    cdef uint32_t min_score = 1 if kmer_length == 3 else 0

    while i < queries._chains.size():
        groups += 1

        group_length = 0
        scores_length = 0

        for j in range(i, queries._chains.size()):

            length = queries._chains[j].get().length() + max_target_length - 2 * kmer_length + 1

            if scores_length + length > max_scores_length and group_length > 0:
                break

            scores_length += length
            group_length += 1

        hash_ = sword.hash.createHash(queries._chains, i, group_length, kmers._kmers)

        for j in range(database_start, database_end):

            for k in range(group_length):
                score_lengths[k] = queries._chains[i + k].get().length() + database._chains[j].get().length() - 2 * kmer_length + 1
                score_starts[k + 1] = score_starts[k] + score_lengths[k]

            sequence = database._chains[j].get().data()
            kmer = sequence[0]

            for k in range(1, kmer_offset):
                kmer = (kmer << kProtBits) | sequence[k]

            max_diag_id = database._chains[j].get().length() - kmer_length
            for k in range(kmer_offset, sequence.size()):
                kmer = ((kmer << kProtBits) | sequence[k]) & del_mask
                hash_.get().hits(begin, end, kmer)

                while begin != end:
                    diagonal = max_diag_id + kmer_offset - k + dereference(begin).position() + score_starts[dereference(begin).id()]
                    scores[diagonal] += 1
                    if max_score[dereference(begin).id()] < scores[diagonal]:
                        max_score[dereference(begin).id()] = scores[diagonal]
                    preincrement(begin)

            for k in range(group_length):
                if max_score[k] <= min_score:
                    continue

                id_ = queries._chains[i + k].get().id()
                flag = entries_part[id_].size() < max_candidates and entries_found[k] < max_candidates

                if flag or max_score[k] >= min_entry_score[k]:
                    entries_part[id_].emplace_back(_ChainEntry(database._chains[j].get().id(), max_score[k]))
                    if min_entry_score[k] > max_score[k]:
                        min_entry_score[k] = max_score[k]

            for k in range(group_length):
                if max_score[k] == 0:
                    continue
                
                max_score[k] = 0
                libcpp.algorithm.fill_n(&scores[0] + score_starts[k], score_lengths[k], 0)

        for k in range(group_length):
            id_ = queries._chains[i + k].get().id()
            dst[id_].insert( dst[id_].end(), entries_part[id_].begin(), entries_part[id_].end() )
            entries_part[id_].clear()

            libcpp.algorithm.stable_sort(dst[id_].begin(), dst[id_].end(), chainEntryDataKey);

            if dst[id_].size() > max_candidates:
                dst[id_].resize(max_candidates)

        i += group_length

    return dst

cpdef object search_database(
    str database_path,
    str queries_path,
    uint32_t kmer_length = 3,
    uint32_t max_candidates = 30000,
    ScoreMatrix score_matrix = ScoreMatrix(),
    uint32_t score_threshold = 13,
):
    cdef vector[uint32_t] tasks
    cdef ChainSet         queries
    cdef ChainSet         db      = ChainSet()
    cdef Reader           dreader = Reader(database_path)
    cdef Reader           qreader = Reader(queries_path)

    # read queries
    queries = qreader.read()
    # libcpp.algorithm.sort(queries._chains.begin(), queries._chains.end(), chainLengthKey)

    # create initial kmers
    cdef Kmers kmers = Kmers(score_matrix, kmer_length, score_threshold)

    # create result entries
    cdef _ChainEntrySet entries# = ChainEntrySet(queries._chains.size())

    while True:

        status = sword.chain.createChainSetPart(db._chains, dreader._reader, 0)
        tasks = preprocess_database(db, 1)

        for i in range(tasks.size() - 1):
            entries = score_chains(
                queries,
                db,
                max_candidates,
                tasks[i],
                tasks[i+1],
                kmers,
            )
   
        if not status:
            break

    # return [[(e.chain_idx(), e.data()) for e in entries[i]] for i in range(len(queries))]

    cdef _Indexes dst = _Indexes(queries._chains.size())
    for i in range(queries._chains.size()):
        dst[i].reserve(entries[i].size())

        for it in entries[i]:
            dst[i].emplace_back(it.chain_idx())

        libcpp.algorithm.sort(dst[i].begin(), dst[i].end())

    return list(dst)



