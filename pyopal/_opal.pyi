import typing

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

SEARCH_MODE = Literal["score", "end", "full"]
SEARCH_OVERFLOW = Literal["simple", "buckets"]
SEARCH_ALGORITHM = Literal["nw", "hw", "ov", "sw"]
COVERAGE_REFERENCE = Literal["query", "target"]

class ScoreMatrix:
    @classmethod
    def aa(cls) -> ScoreMatrix: ...
    def __init__(
        self, alphabet: str, matrix: typing.Iterable[typing.Iterable[int]]
    ) -> None: ...
    def __repr__(self) -> str: ...
    def __reduce__(self) -> typing.Tuple[object, ...]: ...
    @property
    def alphabet(self) -> str: ...
    @property
    def matrix(self) -> typing.List[typing.List[int]]: ...

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

class Database(typing.MutableSequence[str]):
    def __init__(
        self,
        sequences: typing.Iterable[typing.Union[str, bytes, bytearray]] = (),
        score_matrix: typing.Optional[ScoreMatrix] = None,
    ) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> str: ...
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
    @typing.overload
    def search(
        self,
        sequence: typing.Union[str, bytes, bytearray],
        *,
        gap_open: int = 3,
        gap_extend: int = 1,
        mode: Literal["end"] = "end",
        overflow: SEARCH_OVERFLOW = "buckets",
        algorithm: SEARCH_ALGORITHM = "sw",
    ) -> typing.Sequence[EndResult]: ...
    @typing.overload
    def search(
        self,
        sequence: typing.Union[str, bytes, bytearray],
        *,
        gap_open: int = 3,
        gap_extend: int = 1,
        mode: Literal["full"] = "full",
        overflow: SEARCH_OVERFLOW = "buckets",
        algorithm: SEARCH_ALGORITHM = "sw",
    ) -> typing.Sequence[FullResult]: ...
    @typing.overload
    def search(
        self,
        sequence: typing.Union[str, bytes, bytearray],
        *,
        gap_open: int = 3,
        gap_extend: int = 1,
        mode: SEARCH_MODE = "score",
        overflow: SEARCH_OVERFLOW = "buckets",
        algorithm: SEARCH_ALGORITHM = "sw",
    ) -> typing.Sequence[ScoreResult]: ...
