import contextlib
import functools
import multiprocessing.pool
import os
import typing

from scoring_matrices import ScoringMatrix

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

from .lib import (
    Aligner,
    BaseDatabase,
    Database,
    ScoreResult,
    FullResult,
    EndResult,
)

ALIGN_MODE = Literal["score", "end", "full"]
ALIGN_OVERFLOW = Literal["simple", "buckets"]
ALIGN_ALGORITHM = Literal["nw", "hw", "ov", "sw"]

T = typing.TypeVar("T")

@contextlib.contextmanager
def nullcontext(enter_result: T) -> typing.Iterator[T]: ...


@typing.overload
def align(
    query: typing.Union[str, bytes, bytearray],
    database: typing.Union[BaseDatabase, typing.Iterable[typing.Union[str, bytes, bytearray]]],
    scoring_matrix: typing.Union[ScoringMatrix, str, None] = None,
    *,
    gap_open: int = 3,
    gap_extend: int = 1,
    mode: Literal["end"],
    overflow: Literal["simple", "buckets"] = "buckets",
    algorithm: Literal["nw", "hw", "ov", "sw"] = "sw",
    threads: int = 0,
) -> typing.Iterator[EndResult]:
    ...

@typing.overload
def align(
    query: typing.Union[str, bytes, bytearray],
    database: typing.Union[BaseDatabase, typing.Iterable[typing.Union[str, bytes, bytearray]]],
    scoring_matrix: typing.Union[ScoringMatrix, str, None] = None,
    *,
    gap_open: int = 3,
    gap_extend: int = 1,
    mode: Literal["full"],
    overflow: Literal["simple", "buckets"] = "buckets",
    algorithm: Literal["nw", "hw", "ov", "sw"] = "sw",
    threads: int = 0,
) -> typing.Iterator[FullResult]:
    ...

@typing.overload
def align(
    query: typing.Union[str, bytes, bytearray],
    database: typing.Union[BaseDatabase, typing.Iterable[typing.Union[str, bytes, bytearray]]],
    scoring_matrix: typing.Union[ScoringMatrix, str, None] = None,
    *,
    gap_open: int = 3,
    gap_extend: int = 1,
    mode: Literal["score", "end", "full"] = "score",
    overflow: Literal["simple", "buckets"] = "buckets",
    algorithm: Literal["nw", "hw", "ov", "sw"] = "sw",
    threads: int = 0,
    pool: typing.Optional[multiprocessing.pool.ThreadPool] = None
) -> typing.Iterator[ScoreResult]:
    ...