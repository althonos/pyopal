# noqa: D104
from ._version import __version__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__all__ = [
    "Alphabet",
    "Aligner",
    "Database",
    "ScoreMatrix",
    "ScoreResult",
    "EndResult",
    "FullResult",
    "align",
]

import functools
import multiprocessing.pool
import os
import typing

from . import lib
from .lib import (
    Alphabet,
    Aligner,
    Database,
    ScoreMatrix,
    ScoreResult,
    EndResult,
    FullResult,
)

if typing.TYPE_CHECKING:

    try:
        from typing import Literal
    except ImportError:
        from typing_extensions import Literal  # type: ignore

    ALIGN_MODE = Literal["score", "end", "full"]
    ALIGN_OVERFLOW = Literal["simple", "buckets"]
    ALIGN_ALGORITHM = Literal["nw", "hw", "ov", "sw"]

# Use the library documentation
__doc__ = lib.__doc__

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pyopal.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )

def align(
    query: typing.Union[str, bytes, bytearray],
    database: typing.Union[Database, typing.Iterable[typing.Union[str, bytes, bytearray]]],
    score_matrix: typing.Optional[ScoreMatrix] = None,
    *,
    gap_open: int = 3,
    gap_extend: int = 1,
    mode: "ALIGN_MODE" = "score",
    overflow: "ALIGN_OVERFLOW" = "buckets",
    algorithm: "ALIGN_ALGORITHM" = "sw",
    threads: int = 0,
) -> typing.Iterator[ScoreResult]:
    """Align the query sequence to every database sequence in parallel.

    Arguments:
        query (`str` or byte-like object): The sequence to query the
            database with.
        database (iterable of `str` or byte-like objects): The database 
            sequences to align the query to.
        score_matrix (`~pyopal.ScoreMatrix`): The scoring matrix
            to use for the alignment.

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
            of threads reported by `os.cpu_count`. If one given, use
            the main threads for aligning, otherwise spawns a 
            `multiprocessing.pool.ThreadPool`.

    Yields:
        `~pyopal.ScoreResult`: Results for the alignment of the query
        to each target sequence in the database. The actual type depends 
        on the requested ``mode``: it will be `ScoreResult` for mode 
        ``score``, `EndResult` for mode ``end`` and `FullResult` for 
        mode ``full``.

    Example:
        >>> targets = ["AACCGCTG", "ATGCGCT", "TTATTACG"]
        >>> for result in pyopal.align("ACCTG", targets, gap_open=2):
        ...     print(result.score, targets[result.target_index])
        41 AACCGCTG
        31 ATGCGCT
        23 TTATTACG

    """
    # derive default parameters
    if threads == 0:
        threads = os.cpu_count() or 1
    if score_matrix is None:
        score_matrix = ScoreMatrix.aa()
    if not isinstance(database, Database):
        database = Database(database, score_matrix.alphabet)

    # avoid using more threads than necessary
    if threads > len(database):
        threads = len(database) or 1

    # align query to database
    aligner = Aligner(score_matrix, gap_open=gap_open, gap_extend=gap_extend)
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
        # cut the database in chunks of similar length
        chunk_length = len(database) // threads
        with multiprocessing.pool.ThreadPool(threads) as pool:
            align = functools.partial(aligner._align_slice, query, database, mode=mode, overflow=overflow, algorithm=algorithm) # type: ignore
            chunk_hits = pool.imap(
                lambda x: align(start=x, end=x + chunk_length),
                range(0, len(database), chunk_length),    
            )
            for hits in chunk_hits:
                yield from hits
