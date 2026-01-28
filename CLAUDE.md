# CLAUDE.md - Project Context for AI Assistants

## Project Overview

ChromDetect is a Python toolkit for genome assembly classification, validation, and quality control. It provides:

1. **Scaffold Classification** - Classify scaffolds based on naming conventions and size heuristics
2. **Assembly Validation** - Validate FASTA files against NCBI assembly reports
3. **Karyotype Checking** - Verify chromosome counts against expected karyotypes
4. **Name Standardization** - Convert between naming conventions (UCSC, Ensembl, RefSeq, GenBank)
5. **Version Tracking** - Compare assembly versions and track scaffold changes
6. **Multi-Assembly QC Dashboard** - Generate comparative quality reports

**This is a utility tool, not a validator.** It uses pattern matching and size heuristics - it cannot detect misassemblies or verify sequence correctness.

## Key Commands

```bash
# Run tests
pytest

# Run tests with coverage
pytest --cov=chromdetect --cov-report=html

# Type checking
mypy chromdetect

# Linting
ruff check chromdetect

# Install in dev mode
pip install -e ".[dev]"
```

## Architecture

- `chromdetect/core.py` - Main classification logic, ScaffoldInfo/AssemblyStats dataclasses
- `chromdetect/patterns.py` - Regex patterns for scaffold name matching
- `chromdetect/cli.py` - Command-line interface
- `chromdetect/compare.py` - Assembly comparison functionality
- `chromdetect/html_report.py` - HTML report generation
- `chromdetect/assembly_report.py` - NCBI assembly report parsing
- `chromdetect/validation.py` - FASTA validation against assembly reports
- `chromdetect/karyotype.py` - Karyotype database and validation
- `chromdetect/standardize.py` - Naming convention conversion
- `chromdetect/version.py` - Assembly version tracking
- `chromdetect/dashboard.py` - Multi-assembly QC dashboard

## Design Decisions

### Removed Features (v0.5.0)
The following features were removed because they were scientifically inadequate:

1. **Telomere detection** - Simple regex counting of TTAGGG repeats is not reliable for real telomere identification. Users should use TRF or RepeatMasker.

2. **Centromere detection** - AT-rich pattern matching produces too many false positives. Real centromere identification requires specialized tools.

3. **Quality score** - Arbitrary weighted heuristic with no validation against ground-truth datasets. Could mislead users.

### New Features (v0.6.0)
1. **Assembly Report Validation** - Validate FASTA against NCBI reports
2. **Karyotype Consistency Checking** - 29 species database
3. **Assembly Standardization** - UCSC/Ensembl/RefSeq/GenBank conversion
4. **Version Tracking** - Detect promotions/demotions/merges/splits
5. **Multi-Assembly Dashboard** - Comparative QC with HTML charts

### Scope Boundaries
ChromDetect intentionally does NOT:
- Validate assembly correctness
- Detect misassemblies or chimeric scaffolds
- Compare against reference genomes
- Perform synteny analysis

## Version Management

Update version in three places:
1. `chromdetect/__init__.py`
2. `pyproject.toml`
3. `CITATION.cff`

## Citation

DOI: 10.5281/zenodo.17945062

## Testing

386 tests covering:
- Pattern matching for 14+ naming conventions
- Size-based classification
- CLI functionality
- Output formats (JSON, TSV, BED, GFF, HTML)
- Assembly comparison
- NCBI report parsing
- Assembly report validation
- Karyotype validation (29 species)
- Name standardization and NCBI compliance
- Version tracking and comparison
- Multi-assembly QC dashboard
