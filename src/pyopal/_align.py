import contextlib
import functools
import multiprocessing.pool
import os

from scoring_matrices import ScoringMatrix

from .lib import (
    Aligner,
    BaseDatabase,
    Database,
    ScoreResult,
    FullResult,
    EndResult,
)

@contextlib.contextmanager
def nullcontext(enter_result):
    """Return a context manager that returns its input and does nothing.

    Adapted from `contextlib.nullcontext` for backwards compatibility 
    with Python 3.6.

    """
    yield enter_result


def align(
    query,
    database,
    scoring_matrix = None,
    *,
    gap_open = 3,
    gap_extend = 1,
    mode = "score",
    overflow = "buckets",
    algorithm = "sw",
    threads = 0,
    pool = None,
    ordered = False,
):
    """Align the query sequence to every database sequence in parallel.

    Arguments:
        query (`str` or byte-like object): The sequence to query the
            database with.
        database (iterable of `str` or byte-like objects): The database 
            sequences to align the query to. 
        scoring_matrix (`~scoring_matrices.ScoringMatrix` or `str`): The 
            scoring matrix to use for the alignment, either as a 
            `ScoringMatrix` object, or as the name of a matrix to load
            with the `ScoringMatrix.from_name` class method.

    Keyword Arguments:
        gap_open(`int`): The gap opening penalty :math:`G` for
            scoring the alignments.
        gap_extend (`int`): The gap extension penalty :math:`E`
            for scoring the alignments.
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
        threads (`int`): The number of threads to use for aligning
            sequences in parallel. If zero is given, uses the number
            of CPUs reported by `os.cpu_count`. If one given, use
            the main threads for aligning, otherwise spawns a 
            `multiprocessing.pool.ThreadPool`.
        pool (`multiprocessing.pool.ThreadPool`): A running pool 
            instance to use for parallelization. Useful for reusing 
            the same pool across several calls of `~pyopal.align`. 
            If `None` given, spawns a new pool based on the ``threads``
            argument.
        ordered (`bool`): Whether the results should be returned in
            the same order as the database sequences. Internally
            switches the code to use `ThreadPool.imap` instead of
            `ThreadPool.imap_unordered`, which can have an impact
            on performance.

    Yields:
        `~pyopal.ScoreResult`: Results for the alignment of the query
        to each target sequence in the database. The actual type depends 
        on the requested ``mode``: it will be `ScoreResult` for mode 
        ``score``, `EndResult` for mode ``end`` and `FullResult` for 
        mode ``full``.

    Hint:
        Consider storing the database sequences into a `~pyopal.Database` 
        object if you are querying the same sequences more than once 
        to avoid the overhead added by sequence encoding.

    Example:
        >>> targets = ["AACCGCTG", "ATGCGCT", "TTATTACG"]
        >>> for res in pyopal.align("ACCTG", targets, gap_open=2, ordered=True):
        ...     print(res.score, targets[res.target_index])
        41 AACCGCTG
        31 ATGCGCT
        23 TTATTACG

    .. versionadded:: 0.5.0

    """
    # derive default parameters
    if threads == 0:
        threads = os.cpu_count() or 1
    if scoring_matrix is None:
        scoring_matrix = Aligner._DEFAULT_SCORING_MATRIX
    elif isinstance(scoring_matrix, str):
        scoring_matrix = ScoringMatrix.from_name(scoring_matrix)
    elif not isinstance(scoring_matrix, ScoringMatrix):
        ty = type(scoring_matrix).__name__
        raise TypeError(f"expected str or ScoringMatrix, got {ty}")
    if not isinstance(database, BaseDatabase):
        database = Database(database, scoring_matrix.alphabet)

    # avoid using more threads than necessary
    if threads > len(database):
        threads = len(database) or 1

    # align query to database
    aligner = Aligner(scoring_matrix, gap_open=gap_open, gap_extend=gap_extend)
    if threads == 1:
        # use main thread if a single thread is needed
        yield from aligner.align(
            query, 
            database, 
            mode=mode, 
            overflow=overflow, 
            algorithm=algorithm
        )
    else:
        # create a pool if none given
        if pool is None:
            pool_context = multiprocessing.pool.ThreadPool(threads)
        else:
            pool_context = nullcontext(pool)
        # cut the database in chunks of similar length
        chunk_length = len(database) // threads
        with pool_context as pool:
            align = functools.partial(
                aligner.align, # type: ignore
                query, 
                database, 
                mode=mode, 
                overflow=overflow, 
                algorithm=algorithm
            )
            if not ordered:
                chunk_hits = pool.imap_unordered(
                    lambda x: align(start=x, end=x + chunk_length),
                    range(0, len(database), chunk_length),    
                )
            else:
                chunk_hits = pool.imap(
                    lambda x: align(start=x, end=x + chunk_length),
                    range(0, len(database), chunk_length),    
                )
            for hits in chunk_hits:
                yield from hits
