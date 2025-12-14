# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2024-12-14

### Added
- **CLI Enhancements**:
  - `--verbose` / `-v` flag for detailed processing information
  - `--min-confidence FLOAT` filter to include only scaffolds above confidence threshold
  - `--min-length BP` filter to include only scaffolds above length threshold
  - `--list-patterns` flag to display all supported naming patterns
  - Stdin support: use `-` to read FASTA from stdin (e.g., `cat file.fa | chromdetect -`)
  - `__main__.py` for running as module (`python -m chromdetect`)
- **Input Validation**:
  - FASTA format validation with helpful error messages
  - Detection of empty files, missing headers, and invalid formats
- **API Additions**:
  - `parse_fasta_from_handle()` function for reading from file handles/streams
- **Testing**:
  - 21 new CLI integration tests
  - Tests for all CLI flags, output formats, and error conditions
  - Stdin input testing
  - Total: 100 tests with 65% code coverage

### Changed
- `-v` now means `--verbose` (was `--version`)
- `-V` now means `--version`
- Distinct exit codes: 65 for data errors, 66 for input file errors
- Version string now sourced from `__version__` (single source of truth)

### Fixed
- Version flag now correctly reads from package `__version__`

## [0.1.0] - 2024-12-14

### Added
- Initial release of ChromDetect
- Name-based detection with 15+ regex patterns for common naming conventions:
  - Explicit chromosome names (`chr1`, `chromosome_X`)
  - VGP-style scaffolds (`Super_scaffold_1`, `SUPER_1`)
  - Linkage groups (`LG1`, `LG_X`)
  - NCBI accessions (`NC_*`, `CM*`)
  - HiC/RaGOO scaffolder outputs
- Size-based detection using N50 and scaffold size heuristics
- Karyotype-informed classification adjustment
- Command-line interface with multiple output formats:
  - Human-readable summary (default)
  - JSON for programmatic use
  - TSV for spreadsheet tools
- Python API with `parse_fasta()` and `classify_scaffolds()` functions
- Support for gzipped FASTA files
- Assembly statistics (N50, N90, GC content, total length)
- Comprehensive test suite
- Examples for basic usage, batch processing, and karyotype detection

[Unreleased]: https://github.com/shandley/chromdetect/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/shandley/chromdetect/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/shandley/chromdetect/releases/tag/v0.1.0
