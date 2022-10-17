# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyopal/compare/v0.2.0...HEAD


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
