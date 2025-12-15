# CLAUDE.md - Project Context for AI Assistants

## Project Overview

ChromDetect is a Python utility that classifies scaffolds in genome assemblies as chromosomes, unlocalized, or unplaced sequences based on naming conventions and size heuristics.

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
- `chromdetect/fasta.py` - FASTA parsing
- `chromdetect/assembly_report.py` - NCBI assembly report parsing

## Design Decisions

### Removed Features (v0.5.0)
The following features were removed because they were scientifically inadequate:

1. **Telomere detection** - Simple regex counting of TTAGGG repeats is not reliable for real telomere identification. Users should use TRF or RepeatMasker.

2. **Centromere detection** - AT-rich pattern matching produces too many false positives. Real centromere identification requires specialized tools.

3. **Quality score** - Arbitrary weighted heuristic with no validation against ground-truth datasets. Could mislead users.

### Scope Boundaries
ChromDetect intentionally does NOT:
- Validate assembly correctness
- Detect misassemblies or chimeric scaffolds
- Compare against reference genomes
- Perform synteny analysis

It simply answers: "Which scaffolds in this FASTA are likely chromosomes based on their names and sizes?"

## Version Management

Update version in three places:
1. `chromdetect/__init__.py`
2. `pyproject.toml`
3. `CITATION.cff`

## Citation

DOI: 10.5281/zenodo.17945062

## Testing

201 tests covering:
- Pattern matching for 14+ naming conventions
- Size-based classification
- CLI functionality
- Output formats (JSON, TSV, BED, GFF, HTML)
- Assembly comparison
- NCBI report parsing
