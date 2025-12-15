# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] - 2024-12-15

### Removed
- **Telomere Detection**: Removed `--detect-telomeres` flag and related functionality. The simple repeat pattern matching was not scientifically rigorous enough for reliable telomere identification. Users needing telomere analysis should use dedicated tools (e.g., TRF, RepeatMasker).
- **Centromere Detection**: Removed `--detect-centromeres` flag and related functionality. The simplified alpha-satellite pattern matching was inadequate for real centromere identification.
- **Quality Score**: Removed the heuristic quality score from AssemblyStats. The score was not validated against ground-truth datasets and could mislead users about assembly quality.

### Changed
- ScaffoldInfo no longer includes telomere/centromere fields
- AssemblyStats no longer includes telomere_count, t2t_count, centromere_count, or quality_score
- HTML reports simplified to show classification statistics only
- Assembly comparison no longer reports telomere/centromere/quality metrics
- classify_fasta() and compare_fasta_files() simplified API

### Testing
- Reduced to 201 tests after removing telomere/centromere/quality score tests

## [0.4.0] - 2024-12-14

### Added
- **Custom Pattern Support**:
  - `--patterns FILE` flag to load custom patterns from YAML or JSON
  - Merge custom patterns with built-in patterns (prepend for priority)
  - Support for custom chromosome, unlocalized, and fragment patterns

- **NCBI Assembly Report Integration**:
  - `--assembly-report FILE` flag for authoritative scaffold classification
  - Parse NCBI assembly_report.txt files for GenBank/RefSeq assemblies
  - Automatic chromosome count from NCBI reports

- **Assembly Comparison Mode**:
  - `--compare FASTA2` flag to compare two assemblies side-by-side
  - Statistics comparison, chromosome overlap analysis
  - Size differences and classification changes detection
  - Summary and TSV output formats for comparisons

- **HTML Report Output**:
  - `--format html` for visual reports with embedded SVG charts
  - Pie chart for scaffold classification distribution
  - Bar chart for top scaffold sizes
  - Scaffold table with classification details
  - Fully standalone HTML (no external dependencies)

- **Convenience API Functions**:
  - `classify_fasta()` - One-step classification from FASTA file path
  - `compare_fasta_files()` - One-step comparison of two FASTA files
  - Simplifies common workflows without needing to call parse_fasta() separately

- **Per-Scaffold GC Content**:
  - GC content calculated for each scaffold (not just assembly-wide)
  - `gc_content` field added to ScaffoldInfo
  - GC content shown in summary, TSV, JSON, and HTML outputs

### Changed
- classify_scaffolds() now accepts `custom_patterns` and `assembly_report` parameters
- Summary output now shows GC content per scaffold in top 20 list

## [0.3.0] - 2024-12-14

### Added
- **New Output Formats**:
  - BED format (`--format bed`) for integration with bedtools, IGV, and genomics pipelines
  - GFF3 format (`--format gff`) with classification metadata in attributes
- **Chromosome Sequence Extraction**:
  - `--extract-chromosomes FILE` flag to extract chromosome sequences to a new FASTA file
  - `write_fasta()` function in Python API for writing FASTA files
  - `parse_fasta(keep_full_sequence=True)` option to retain full sequences
- **Batch Processing**:
  - `--batch DIR` flag to process all FASTA files in a directory
  - Automatic output file generation per assembly
  - Batch summary TSV with overview of all processed assemblies
- **API Additions**:
  - `format_bed()` function for BED format output
  - `format_gff()` function for GFF3 format output
  - `write_fasta()` function for FASTA file writing

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
- Name-based detection with 14 regex patterns for common naming conventions:
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

[Unreleased]: https://github.com/shandley/chromdetect/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/shandley/chromdetect/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/shandley/chromdetect/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/shandley/chromdetect/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/shandley/chromdetect/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/shandley/chromdetect/releases/tag/v0.1.0
