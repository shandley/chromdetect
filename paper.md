---
title: 'ChromDetect: A utility for classifying scaffolds in genome assemblies'
tags:
  - Python
  - genomics
  - bioinformatics
  - genome assembly
authors:
  - name: Scott A. Handley
    orcid: 0000-0002-2143-6570
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Pathology and Immunology, Washington University School of Medicine, St. Louis, Missouri, USA
    index: 1
date: 14 December 2024
bibliography: paper.bib
---

# Summary

Genome assemblies lack standardized naming conventions for scaffolds. Different assemblers and scaffolding tools produce vastly different nomenclature—`chr1`, `Super_scaffold_1`, `LG_1`, `HiC_scaffold_12`, `NC_000001.11`—making it difficult to programmatically identify which scaffolds represent chromosomes.

ChromDetect is a Python utility that classifies scaffolds in genome assemblies as chromosomes, unlocalized, or unplaced sequences. It uses pattern matching against common naming conventions and size-based heuristics to automate a task that would otherwise require manual inspection or custom scripts.

# Statement of Need

When working with genome assemblies, researchers often need to:

- Filter an assembly to extract only chromosome-level scaffolds
- Count how many chromosomes an assembly contains
- Compare scaffold classifications across different assembly versions
- Generate chromosome lists for downstream tools

These tasks are straightforward when scaffolds are named `chr1`, `chr2`, etc., but become tedious when assemblies use conventions like `Super_scaffold_1`, `HiC_scaffold_12`, or `scaffold_28_cov50`.

ChromDetect automates this classification. It is not an assembly validator—tools like QUAST [@gurevich2013quast] and Merqury handle quality assessment. ChromDetect simply answers: "Which scaffolds in this FASTA file are chromosomes?"

# Approach

ChromDetect uses simple, interpretable rules:

1. **Name matching**: Regular expressions match 14 common naming patterns (`chr1`, `Super_scaffold_1`, `LG_1`, `NC_*`, `HiC_scaffold_*`, etc.)

2. **Size heuristics**: Scaffolds above a threshold (default 10 Mb) are likely chromosomes

3. **Karyotype adjustment** (optional): If you know the expected chromosome count, ChromDetect adjusts classifications accordingly

4. **NCBI report parsing** (optional): Uses official NCBI assembly reports for authoritative classification

# Features

ChromDetect provides:

- Command-line interface and Python API
- Multiple output formats (JSON, TSV, BED, GFF3, HTML)
- Assembly comparison mode for evaluating assembly improvements
- Batch processing of multiple assemblies
- Zero external dependencies (pure Python, works with Python 3.9-3.12)
- Comprehensive test suite (201 tests)

# Example Usage

```bash
# Basic analysis
chromdetect assembly.fasta

# Generate visual HTML report
chromdetect assembly.fasta --format html -o report.html

# Compare two assemblies
chromdetect assembly_v1.fasta --compare assembly_v2.fasta

# Use known karyotype for adjustment
chromdetect assembly.fasta --karyotype 24
```

# Implementation

ChromDetect is implemented in Python with no external dependencies. The core classification algorithm combines confidence scores from multiple detection methods using priority rules: strong name matches take precedence, followed by size-based detection for ambiguous cases.

The tool is available on PyPI (`pip install chromdetect`) and GitHub. Continuous integration ensures compatibility across Python 3.9-3.12 on Linux, macOS, and Windows platforms.

# Limitations

ChromDetect is designed for scaffold classification, not assembly validation. Key limitations include:

- **Pattern-dependent**: Classification relies on naming conventions; unusual naming schemes require custom patterns
- **No misassembly detection**: Cannot identify structural errors, chimeric scaffolds, or sequence inaccuracies
- **No reference comparison**: Does not compare against reference genomes or identify missing chromosomes

ChromDetect should complement, not replace, comprehensive assembly QC tools like QUAST [@gurevich2013quast] or Merqury for validation workflows.

# Acknowledgments

We thank the Vertebrate Genomes Project, T2T Consortium, and genome assembly community for discussions on scaffold naming conventions and classification challenges.

# References
