# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyopal/compare/v0.5.0-a2...HEAD


## [v0.5.0-a2] - 2024-01-12
[v0.5.0-a2]: https://github.com/althonos/pyopal/compare/v0.5.0-a1...v0.5.0-a2

### Added
- `pyopal.Alphabet` class for handling sequence encoding outside of `ScoreMatrix`.

### Changed
- Make internal extension public as `pyopal.lib` submodule.


## [v0.5.0-a1] - 2024-01-12
[v0.5.0-a1]: https://github.com/althonos/pyopal/compare/v0.4.2...v0.5.0-a1

### Changed
- Store `Database` sequences using `shared_ptr` for fast slicing and extraction of a database subset.
- Use templates to generate platform-specific Cython files.

### Fixed
- SSE2 not being detected on x86-64, where it is mandatory.
- Mark Opal search functions as `noexcept` in Cython to avoid exception checks.

### Added
- Support for compiling and running package on Windows.
- Support for Python 3.12.


## [v0.4.2] - 2023-08-29
[v0.4.2]: https://github.com/althonos/pyopal/compare/v0.4.1...v0.4.2

### Fixed
- Segmentation fault caused by out-of-bounds access in Opal ([#5](https://github.com/althonos/pyopal/issues/5)).

### Changed
- Use Cython `binding` directive to improve generated documentation. 


## [v0.4.1] - 2023-08-29
[v0.4.1]: https://github.com/althonos/pyopal/compare/v0.4.0...v0.4.1

### Fixed
- `FullResult.__repr__` not returning a roundtripping string ([#4](https://github.com/althonos/pyopal/pull/4), by [@valentynbez](https://github.com/valentynbez)).
- `Database.search` overflowing for long sequences with non-`sw` algorithms ([#3](https://github.com/althonos/pyopal/issues/3)).

### Changed
- Make `Database.search` raise an `OverflowError` instead of a `RuntimeError` on score overflow.


## [v0.4.0] - 2023-07-21
[v0.4.0]: https://github.com/althonos/pyopal/compare/v0.3.0...v0.4.0

### Changed
- Bumped Cython dependency to `v3.0`.

### Fixed
- PyPy 3.9 builds failing on missing `PyInterpreterState_GetID`.


## [v0.3.0] - 2023-05-30
[v0.3.0]: https://github.com/althonos/pyopal/compare/v0.2.0...v0.3.0

### Added
- `Database.extract` method to subset a database given a sequence of indices.
- `Database.mask` method to subset a database given a boolean mask.

### Changed
- Replaced `cpu-features` library with `archspec` Python package for runtime detection of CPU features.

### Fixed
- Segmentation fault in alignment reconstruction code for Needleman-Wunsch algorithm ([#1](https://github.com/althonos/pyopal/issues/1)).
- Erroneous error message in `Database.search` on invalid `overflow` value ([#2](https://github.com/althonos/pyopal/issues/2)).


## [v0.2.0] - 2022-10-17
[v0.2.0]: https://github.com/althonos/pyopal/compare/v0.1.1...v0.2.0

### Added
- `query_length` and `target_length` properties to `FullResult` to the lenghts of the complete query and target sequences.
- `FullResult.coverage` method to compute the coverage of the alignment using either the query or the target as the reference.

### Changed
- Compile Cython extension with `binding=False` to fix rendering of documentation.

### Fixed
- Insertion & deletion symbols being inverted in `FullResult.cigar` strings.


## [v0.1.1] - 2022-10-07
[v0.1.1]: https://github.com/althonos/pyopal/compare/v0.1.0...v0.1.1

### Added
- Buffer protocol implementation to `pyopal.ScoreMatrix`.
- Sphinx documentation hosted on [ReadTheDocs](https://pyopal.readthedocs.io).

### Fixed
- Docstring now showing in the main `pyopal` module.
- `Database.insert` potentially crashing when given negative indexes.
- `Database.search` not listing the right return type.


## [v0.1.0] - 2022-10-06
[v0.1.0]: https://github.com/althonos/pyopal/compare/347b6d3...v0.1.0

Initial release.
