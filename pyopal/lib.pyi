import sys
import typing

from scoring_matrices import ScoringMatrix

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

ALIGN_MODE = Literal["score", "end", "full"]
ALIGN_OVERFLOW = Literal["simple", "buckets"]
ALIGN_ALGORITHM = Literal["nw", "hw", "ov", "sw"]
COVERAGE_REFERENCE = Literal["query", "target"]

# --- Parameters ---------------------------------------------------------------

class Alphabet(typing.Sequence[str]):
    def __init__(self, letters: str = "ARNDCQEGHILKMFPSTWYVBZX*"): ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> str: ...
    def __contains__(self, value: object) -> bool: ...
    def __reduce__(self) -> typing.Tuple[typing.Type[Alphabet], typing.Tuple[str]]: ...
    def __repr__(self) -> str: ...
    def __eq__(self, obj: object) -> bool: ...
    def __str__(self) -> str: ...
    def encode(self, sequence: typing.Union[bytes, bytearray, memoryview]) -> bytes: ...
    def decode(self, encoded: typing.Union[bytes, bytearray, memoryview]) -> str: ...


# --- Sequence storage ---------------------------------------------------------

class BaseDatabase(typing.Sequence[str]):
    @property
    def alphabet(self) -> Alphabet: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> str: ...

class Database(typing.MutableSequence[str], BaseDatabase):
    def __init__(
        self,
        sequences: typing.Iterable[typing.Union[str, bytes, bytearray]] = (),
        alphabet: typing.Union[Alphabet, str, None] = None,
    ) -> None: ...
    @typing.overload
    def __getitem__(self, index: int) -> str: ...
    @typing.overload
    def __getitem__(self, index: slice) -> Database: ...
    @typing.overload
    def __getitem__(self, index: typing.Union[int, slice]) -> typing.Union[str, Database]: ...
    def __setitem__(
        self, index: int, sequence: typing.Union[str, bytes, bytearray]
    ) -> None: ...
    def __delitem__(self, index: int) -> None: ...
    def clear(self) -> None: ...
    def extend(
        self, sequences: typing.Iterable[typing.Union[str, bytes, bytearray]]
    ) -> None: ...
    def append(self, sequence: typing.Union[str, bytes, bytearray]) -> None: ...
    def reverse(self) -> None: ...
    def insert(
        self, index: int, sequence: typing.Union[str, bytes, bytearray]
    ) -> None: ...
    def mask(self, bitmask: typing.Sequence[bool]) -> Database: ...
    def extract(self, indices: typing.Sequence[int]) -> Database: ...

# --- Aligner ------------------------------------------------------------------

class ScoreResult:
    def __init__(self, target_index: int, score: int) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def target_index(self) -> int: ...
    @property
    def score(self) -> int: ...

class EndResult(ScoreResult):
    def __init__(
        self, target_index: int, score: int, query_end: int, target_end: int
    ) -> None: ...
    @property
    def query_end(self) -> int: ...
    @property
    def target_end(self) -> int: ...

class FullResult(EndResult):
    def __init__(
        self,
        target_index: int,
        score: int,
        query_end: int,
        target_end: int,
        query_start: int,
        target_start: int,
        query_length: int,
        target_length: int,
        alignment: str,
    ) -> None: ...
    @property
    def query_start(self) -> int: ...
    @property
    def target_start(self) -> int: ...
    @property
    def query_length(self) -> int: ...
    @property
    def target_length(self) -> int: ...
    @property
    def alignment(self) -> str: ...
    def cigar(self) -> str: ...
    def identity(self) -> float: ...
    def coverage(self, reference: COVERAGE_REFERENCE = "query") -> float: ...

class Aligner:
    def __init__(
        self,
        scoring_matrix: typing.Union[ScoringMatrix, str, None] = None,
        *,
        gap_open: int = 3,
        gap_extend: int = 10
    ) -> None: ...
    @property
    def scoring_matrix(self) -> ScoringMatrix: ...
    @property
    def alphabet(self) -> Alphabet: ...
    @property
    def gap_open(self) -> int: ...
    @property
    def gap_extend(self) -> int: ...
    @typing.overload
    def align(
        self,
        query: typing.Union[str, bytes, bytearray],
        database: BaseDatabase,
        *,
        mode: Literal["full"] = "full",
        overflow: ALIGN_OVERFLOW = "buckets",
        algorithm: ALIGN_ALGORITHM = "sw",
        start: int = 0,
        end: int = sys.maxsize,
    ) -> typing.Sequence[FullResult]: ...
    @typing.overload
    def align(
        self,
        query: typing.Union[str, bytes, bytearray],
        database: BaseDatabase,
        *,
        mode: Literal["end"] = "end",
        overflow: ALIGN_OVERFLOW = "buckets",
        algorithm: ALIGN_ALGORITHM = "sw",
        start: int = 0,
        end: int = sys.maxsize,
    ) -> typing.Sequence[EndResult]: ...
    @typing.overload
    def align(
        self,
        query: typing.Union[str, bytes, bytearray],
        database: BaseDatabase,
        *,
        mode: ALIGN_MODE = "score",
        overflow: ALIGN_OVERFLOW = "buckets",
        algorithm: ALIGN_ALGORITHM = "sw",
        start: int = 0,
        end: int = sys.maxsize,
    ) -> typing.Sequence[ScoreResult]: ...
