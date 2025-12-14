# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/shandley/chromdetect/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/shandley/chromdetect/releases/tag/v0.1.0
