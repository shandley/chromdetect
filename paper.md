---
title: 'ChromDetect: Automated scaffold classification and validation for genome assemblies'
tags:
  - Python
  - genomics
  - bioinformatics
  - genome assembly
  - scaffold classification
authors:
  - name: Scott A. Handley
    orcid: 0000-0002-2143-6570
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Pathology and Immunology, Washington University School of Medicine, St. Louis, Missouri, USA
    index: 1
date: 28 January 2026
bibliography: paper.bib
---

# Abstract

**Summary:** ChromDetect is a Python toolkit for classifying scaffolds in genome assemblies as chromosomes, unlocalized, or unplaced sequences. Unlike existing assembly statistics tools that compute metrics like N50 without distinguishing scaffold types, ChromDetect uses pattern matching against 14 common naming conventions and size-based heuristics to automate classification. The toolkit provides assembly validation against NCBI reports, karyotype verification for 29 species, name standardization across four conventions (UCSC, Ensembl, RefSeq, GenBank), version tracking between assemblies, and multi-assembly quality control dashboards.

**Availability and Implementation:** ChromDetect is implemented in Python with zero external dependencies, supporting Python 3.9-3.12 on Linux, macOS, and Windows. Source code is available at https://github.com/shandley/chromdetect under the MIT license. The package is installable via PyPI (`pip install chromdetect`).

**Contact:** handley.scott@gmail.com

**Supplementary Information:** Documentation and examples available at https://github.com/shandley/chromdetect

# Introduction

Genome assemblies use inconsistent naming conventions for scaffolds. Assembly tools produce varied nomenclature—`chr1`, `Super_scaffold_1`, `LG_1`, `HiC_scaffold_12`, `NC_000001.11`, `scaffold_28_cov50`—reflecting different software pipelines, institutions, and historical practices [@rhie2021towards]. This heterogeneity creates practical challenges when researchers need to identify chromosome-level scaffolds, validate assemblies against expected karyotypes, or standardize names for cross-study comparisons.

Existing tools address assembly quality assessment but not scaffold classification. QUAST [@gurevich2013quast] evaluates assemblies using metrics like N50, misassembly counts, and gene completeness but does not classify individual scaffolds by type. Similarly, gfastats [@formenti2022gfastats] computes comprehensive statistics for FASTA/GFA files but treats all scaffolds uniformly without distinguishing chromosomes from fragments. Assembly-stats and related utilities likewise report aggregate metrics without scaffold-level classification.

ChromDetect fills this gap by providing automated scaffold classification and validation. It answers the question that existing tools do not: "Which scaffolds in this FASTA file represent chromosomes, and do they match expectations?"

# Implementation

ChromDetect is implemented in pure Python with no external dependencies, ensuring straightforward installation and broad compatibility. The toolkit provides both a command-line interface and a Python API.

## Scaffold Classification

Classification uses two complementary approaches. First, regular expressions match scaffolds against 14 common naming patterns including chromosome prefixes (`chr1`, `Chr_1`, `chromosome_1`), super scaffolds (`Super_scaffold_1`, `SUPER_1`), linkage groups (`LG1`, `LG_1`), NCBI accessions (`NC_000001.11`, `CM000663.2`), and assembly tool outputs (`HiC_scaffold_1`, `Scaffold_1_RaGOO`). Second, size-based heuristics classify large scaffolds (default threshold: 10 Mb) as likely chromosomes when name patterns are ambiguous.

Each scaffold receives a classification (chromosome, unlocalized, or unplaced) and confidence score based on match strength. The algorithm prioritizes strong name matches over size-based detection, with user-configurable thresholds.

## Assembly Validation

ChromDetect validates FASTA files against NCBI assembly reports, detecting mismatches in sequence counts, scaffold lengths, and naming. This ensures downloaded or processed assemblies match official specifications. A strict mode enforces zero-tolerance validation for production pipelines.

## Karyotype Verification

A built-in database covering 29 species enables karyotype checking. Users can verify that assemblies contain expected chromosome counts for organisms including human (23 pairs), mouse (20), zebrafish (25), *Drosophila* (4), *Arabidopsis* (5), and *S. cerevisiae* (16). The database includes mammals, fish, invertebrates, plants, and microorganisms commonly used in genomics research.

## Name Standardization

ChromDetect converts scaffold names between four conventions: UCSC (`chr1`, `chrX`), Ensembl (`1`, `X`), RefSeq (`NC_000001.11`), and GenBank (`CM000663.2`). This facilitates data integration across resources using different naming schemes. An NCBI compliance checker validates names against submission requirements.

## Version Tracking

Assembly comparison functionality tracks changes between versions, detecting scaffold promotions (unplaced to chromosome), demotions, splits, and merges. Summary statistics include N50 changes and scaffold count differences, useful for evaluating assembly improvements.

## Quality Control Dashboard

For multi-assembly projects, ChromDetect generates HTML dashboards with interactive visualizations comparing scaffold counts, N50 values, total sizes, and classification distributions across assemblies. Automatic QC flags highlight potential issues such as fragmentation or low chromosome counts.

# Usage Examples

Basic classification:
```bash
chromdetect assembly.fasta
chromdetect assembly.fasta --format json -o results.json
```

Validation and karyotype checking:
```bash
chromdetect assembly.fasta --assembly-report report.txt --validate
chromdetect assembly.fasta --check-karyotype human
```

Name standardization and version comparison:
```bash
chromdetect assembly.fasta --rename ucsc -o standardized.fasta
chromdetect v1.fasta --compare-versions v2.fasta
```

Multi-assembly dashboard:
```bash
chromdetect --dashboard *.fasta -o dashboard.html --format html
```

Python API:
```python
from chromdetect import classify_fasta

results, stats = classify_fasta("assembly.fasta")
chromosomes = [r for r in results if r.classification == "chromosome"]
print(f"Found {len(chromosomes)} chromosomes, N50: {stats.n50/1e6:.1f} Mb")
```

# Comparison with Existing Tools

Figure 1 compares ChromDetect features with related assembly tools. While QUAST [@gurevich2013quast], gfastats [@formenti2022gfastats], and assembly-stats provide assembly statistics, ChromDetect uniquely offers scaffold classification by naming convention, karyotype verification, name standardization, and NCBI report validation. These tools are complementary: QUAST provides reference-based quality assessment including misassembly detection, while gfastats offers GFA format support and sequence manipulation. ChromDetect addresses the distinct problem of scaffold classification and validation.

![**Figure 1. ChromDetect features and validation.** (A) Feature comparison with existing assembly tools. White circles indicate supported features. ChromDetect (highlighted) uniquely provides scaffold classification, karyotype verification, and name standardization. (B) Validation accuracy across eight diverse genome assemblies spanning vertebrates (blue), plants (green), bacteria (purple), and telomere-to-telomere assemblies (red). ChromDetect achieved 100% classification accuracy across 283 scaffolds. Numbers above bars indicate accuracy and scaffold count.](figure1.png)

# Validation

ChromDetect was validated against NCBI assembly reports spanning eight diverse genome assemblies: vertebrates (zebra finch VGP, chicken, human GRCh38, zebrafish), plants (*Arabidopsis thaliana* TAIR10), bacteria (*E. coli* K-12), and the telomere-to-telomere human assembly (T2T-CHM13v2.0). Ground truth classifications were derived from official NCBI assembly reports.

Classification accuracy was **100%** across 283 test scaffolds with diverse naming conventions: SUPER_* (VGP), NC_*/NW_* (RefSeq accessions), numeric chromosomes, accession-based names, and ambiguous contig names (CTG*, NA*). The validation demonstrates accurate detection across vertebrate, plant, and bacterial genomes, as well as both draft (GRCh38) and complete (T2T) human assemblies. The test suite includes 404 unit tests covering pattern matching, size classification, CLI functionality, and all output formats.

# Limitations

ChromDetect relies on naming conventions and size heuristics. It cannot detect misassemblies, sequence errors, or chimeric scaffolds—problems requiring reference-based tools like QUAST or read-mapping approaches like Inspector. Unusual or custom naming schemes may require user-defined patterns via YAML configuration. ChromDetect is designed for classification and validation, not comprehensive assembly QC.

# Conclusion

ChromDetect provides automated scaffold classification, validation, and standardization for genome assemblies. By addressing a gap between assembly statistics tools and comprehensive QC pipelines, it simplifies workflows requiring scaffold-level information. The zero-dependency design and cross-platform support ensure accessibility for diverse research environments.

# Funding

This work was not supported by external funding.

# References
