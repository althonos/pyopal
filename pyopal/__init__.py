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
import os
import multiprocessing.pool

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
    query,
    database,
    score_matrix = None,
    *,
    gap_open = 3,
    gap_extend = 1,
    mode = "score",
    overflow = "buckets",
    algorithm = "sw",
    threads = 0,
):
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
            align = functools.partial(aligner._align_slice, query, database, mode=mode, overflow=overflow, algorithm=algorithm)
            chunk_hits = pool.imap(
                lambda x: align(start=x, end=x + chunk_length),
                range(0, len(database), chunk_length),    
            )
            for hits in chunk_hits:
                yield from hits
        

        



