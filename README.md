<p align="center">
  <img src="chromdetect_hex_logo.jpeg" alt="ChromDetect Logo" width="200"/>
</p>

<h1 align="center">ChromDetect</h1>

<p align="center">
  <a href="https://pypi.org/project/chromdetect/"><img src="https://img.shields.io/pypi/v/chromdetect.svg" alt="PyPI version"></a>
  <a href="https://pypi.org/project/chromdetect/"><img src="https://img.shields.io/pypi/pyversions/chromdetect.svg" alt="Python versions"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://github.com/shandley/chromdetect/actions/workflows/test.yml"><img src="https://github.com/shandley/chromdetect/actions/workflows/test.yml/badge.svg" alt="Tests"></a>
  <a href="https://doi.org/10.5281/zenodo.17945062"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17945062.svg" alt="DOI"></a>
</p>

<p align="center">
  <strong>A toolkit for genome assembly classification, validation, and quality control.</strong>
</p>

---

## Overview

ChromDetect helps you work with genome assemblies by providing six key capabilities:

| Feature | Description |
|---------|-------------|
| **Scaffold Classification** | Identify chromosomes vs unplaced scaffolds based on naming patterns and size |
| **Assembly Validation** | Validate FASTA files against NCBI assembly reports |
| **Karyotype Checking** | Verify chromosome counts against 29 species databases |
| **Name Standardization** | Convert between UCSC, Ensembl, RefSeq, and GenBank conventions |
| **Version Tracking** | Compare assembly versions and detect scaffold changes |
| **QC Dashboard** | Generate comparative reports across multiple assemblies |

## Installation

```bash
pip install chromdetect
```

## Quick Examples

```bash
# Classify scaffolds in an assembly
chromdetect assembly.fasta

# Validate against NCBI report
chromdetect assembly.fasta --assembly-report report.txt --validate

# Check chromosome count for human
chromdetect assembly.fasta --check-karyotype human

# Convert to UCSC naming (chr1, chr2, chrX)
chromdetect assembly.fasta --rename ucsc -o renamed.fasta

# Compare two assembly versions
chromdetect v1.fasta --compare-versions v2.fasta

# Generate QC dashboard for multiple assemblies
chromdetect --dashboard *.fasta -o dashboard.html --format html
```

## Use Cases

### Preparing assemblies for submission

Before submitting to NCBI, check compliance and standardize names:

```bash
# Check if names meet NCBI requirements
chromdetect assembly.fasta --check-compliance

# Rename to standard convention
chromdetect assembly.fasta --rename refseq -o submission_ready.fasta
```

### Quality control across projects

Compare multiple assemblies from different sources:

```bash
# Generate comparative dashboard
chromdetect --dashboard sample1.fa sample2.fa sample3.fa -o qc_report.html --format html
```

### Validating downloaded assemblies

Verify a FASTA matches its NCBI assembly report:

```bash
chromdetect GRCh38.fasta --assembly-report GRCh38_report.txt --validate --strict
```

### Tracking assembly improvements

See what changed between versions:

```bash
chromdetect old_assembly.fasta --compare-versions new_assembly.fasta
```

Output shows promotions, demotions, and metric changes:
```
SCAFFOLD CHANGES:
  Promoted:    2 scaffolds (unplaced → chromosome)
  Unchanged:   1,150 scaffolds
  N50 change:  +6.7 Mb (+14.6%)
```

### Checking species-specific karyotype

Verify your assembly has the expected chromosomes:

```bash
# List available species
chromdetect --list-species

# Check against expected karyotype
chromdetect mouse_assembly.fasta --check-karyotype mouse
```

## Output Formats

| Format | Flag | Use Case |
|--------|------|----------|
| Summary | `--format summary` | Quick terminal inspection (default) |
| JSON | `--format json` | Programmatic processing |
| TSV | `--format tsv` | Spreadsheet analysis |
| HTML | `--format html` | Visual reports with charts |
| BED | `--format bed` | Genomics pipelines (bedtools, etc.) |
| GFF | `--format gff` | Genome browsers |

## Python API

```python
from chromdetect import classify_fasta

# Classify an assembly
results, stats = classify_fasta("assembly.fasta")
print(f"Chromosomes: {stats.chromosome_count}")
print(f"N50: {stats.n50 / 1e6:.1f} Mb")

# Filter to just chromosomes
chromosomes = [r for r in results if r.classification == "chromosome"]
for c in chromosomes:
    print(f"  {c.name}: {c.length:,} bp")
```

Additional modules for specific tasks:

```python
# Validation
from chromdetect.validation import validate_fasta_against_report

# Karyotype checking
from chromdetect.karyotype import validate_karyotype, KaryotypeDatabase

# Name standardization
from chromdetect.standardize import standardize_fasta, check_ncbi_compliance

# Version comparison
from chromdetect.version import compare_fasta_files

# Multi-assembly dashboard
from chromdetect.dashboard import analyze_multiple_assemblies, generate_dashboard_html
```

## Supported Species (Karyotype Database)

ChromDetect includes karyotype data for 29 species:

**Mammals:** Human, mouse, rat, dog, cat, horse, cow, pig, sheep, goat, rabbit, guinea pig

**Other vertebrates:** Chicken, zebrafish, frog

**Invertebrates:** Fruit fly, C. elegans

**Plants:** Arabidopsis, rice, maize, wheat, soybean, tomato

**Microorganisms:** Yeast (S. cerevisiae), E. coli

Use `chromdetect --list-species` to see all available species with chromosome counts.

## Recognized Naming Patterns

ChromDetect automatically recognizes common scaffold naming conventions:

- **Chromosome prefixes:** `chr1`, `Chr_1`, `chromosome_1`, `Chromosome1`
- **Super scaffolds:** `Super_scaffold_1`, `Superscaffold_1`, `SUPER_1`
- **Linkage groups:** `LG1`, `LG_1`, `linkage_group_1`
- **NCBI accessions:** `NC_000001.11`, `CM000663.2`
- **Assembly tools:** `HiC_scaffold_1`, `Scaffold_1_RaGOO`
- **Simple numeric:** `1`, `2`, `X`, `MT`

Custom patterns can be added via YAML configuration files.

## Limitations

ChromDetect uses naming patterns and size heuristics—it cannot:
- Detect misassemblies or sequence errors
- Validate sequence correctness
- Perform synteny or homology analysis

For comprehensive assembly validation, use ChromDetect alongside tools like [QUAST](https://github.com/ablab/quast) or [Merqury](https://github.com/marbl/merqury).

## Citation

If you use ChromDetect in your research, please cite:

```bibtex
@software{chromdetect,
  author = {Handley, Scott A.},
  title = {ChromDetect: A toolkit for genome assembly classification and QC},
  url = {https://github.com/shandley/chromdetect},
  version = {0.6.0},
  doi = {10.5281/zenodo.17945062},
  year = {2025}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
