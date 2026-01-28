# ChromDetect Validation Report

**Date:** January 28, 2026
**ChromDetect Version:** 0.6.0
**Validation Method:** Comparison against NCBI assembly reports

## Overview

This report documents validation of ChromDetect's scaffold classification accuracy
using real genome assemblies from NCBI. Ground truth classifications are derived
from official NCBI assembly reports.

## Test Assemblies

| Assembly | Species | Source | Naming Convention | Scaffolds Tested |
|----------|---------|--------|-------------------|------------------|
| bTaeGut1.4.pri | Zebra finch (*Taeniopygia guttata*) | VGP | SUPER_* | 50 |
| bTaeGut1.4.pri | Zebra finch | VGP | NC_* (RefSeq) | 50 |
| GRCg6a | Chicken (*Gallus gallus*) | GRC | Numeric (1, 2, 3...) | 50 |
| GRCh38.p14 | Human (*Homo sapiens*) | GRC | Numeric + HSCHR* | 50 |
| GRCz11 | Zebrafish (*Danio rerio*) | GRC | Numeric + CTG*/NA* | 50 |

## Results Summary

| Assembly | Naming Convention | Accuracy | Correct | Incorrect |
|----------|-------------------|----------|---------|-----------|
| Zebra finch | SUPER_* (VGP) | **82.0%** | 41 | 9 |
| Zebra finch | NC_* (RefSeq) | **82.0%** | 41 | 9 |
| Chicken | Numeric | **100.0%** | 50 | 0 |
| Human | Numeric | **98.0%** | 49 | 1 |
| Zebrafish | Numeric + ambiguous | **50.0%** | 25 | 25 |
| **OVERALL** | | **82.4%** | 206 | 44 |

## Analysis by Naming Convention

### Standard Naming Conventions (High Accuracy)

**Chicken (100%)** and **Human (98%)** assemblies use standard NCBI/GRC naming:
- Numeric chromosome names: `1`, `2`, `3`, ...
- Clear unlocalized patterns: `chr1_random`, etc.
- ChromDetect correctly identifies these conventions.

### VGP SUPER_* Convention (82%)

Zebra finch VGP assembly uses `SUPER_N` naming for chromosomes:
- ChromDetect correctly identifies SUPER_* scaffolds as chromosomes
- **Issue:** Unplaced scaffolds named `scaffold_N_ctg1` are misclassified
- The `scaffold_*` pattern overlaps with assembly contig naming

### RefSeq NC_* Accessions (82%)

- `NC_*` accessions (assembled chromosomes) correctly identified
- `NW_*` accessions (unplaced scaffolds) misclassified as chromosomes
- ChromDetect lacks specific pattern for NW_* (unplaced) accessions

### Ambiguous Naming (50%)

Zebrafish assembly contains scaffolds with ambiguous names:
- `NA876`, `CTG2805`, `CTG386` - no clear chromosome indicator
- ChromDetect defaults to chromosome classification
- These should be classified as unplaced (unknown)

## Specific Misclassification Patterns

### 1. Unplaced NW_* Accessions
```
NW_024545304.1: predicted=chromosome, expected=unplaced
```
**Recommendation:** Add NW_* pattern to FRAGMENT_PATTERNS

### 2. Assembly Contig Names
```
scaffold_104_ctg1: predicted=chromosome, expected=unplaced
```
**Recommendation:** The `_ctg\d` suffix should trigger unplaced classification

### 3. Ambiguous Short Names
```
NA876: predicted=chromosome, expected=unplaced
CTG2805: predicted=chromosome, expected=unplaced
```
**Recommendation:** Consider "unknown" classification for ambiguous names

## Recommendations for Improvement

1. **Add NW_* pattern to unplaced scaffolds**
   ```python
   FRAGMENT_PATTERNS.append(r'^NW_\d+\.\d+$')
   ```

2. **Strengthen CTG/contig detection**
   - Current: `r'ctg\d*$'`
   - Improvement: `r'(^|_)ctg\d*'` to catch `scaffold_N_ctg1`

3. **Add "unknown" classification**
   - For scaffolds that don't match any pattern and are small
   - Reduces false positive chromosome calls

4. **Consider NCBI accession prefix semantics**
   - NC_* = RefSeq chromosome
   - NW_* = RefSeq unplaced scaffold
   - NT_* = RefSeq genomic contig

## Conclusion

ChromDetect achieves **82.4% overall accuracy** on real genome assemblies.
Accuracy is highest (98-100%) on assemblies using standard NCBI/GRC naming
conventions. Performance degrades with ambiguous or non-standard naming schemes.

The validation identifies specific patterns for improvement, particularly:
- NW_* accession recognition
- Better handling of assembly contig suffixes
- Ambiguous name handling

These findings are documented to guide future development and to set realistic
expectations for users working with diverse assembly types.

## Files

- `validation_results.json` - Detailed JSON results
- `validate_accuracy.py` - Validation script
- `*_report.txt` - NCBI assembly reports
- `*_test.fasta` - Synthetic test FASTA files
