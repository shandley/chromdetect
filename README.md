# ChromDetect

[![PyPI version](https://badge.fury.io/py/chromdetect.svg)](https://badge.fury.io/py/chromdetect)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/verity-project/chromdetect/actions/workflows/test.yml/badge.svg)](https://github.com/verity-project/chromdetect/actions/workflows/test.yml)

**Detect chromosome-level scaffolds in genome assemblies with inconsistent naming conventions.**

## The Problem

Genome assemblies use wildly inconsistent naming conventions for chromosome-level scaffolds:

- `Super_scaffold_1`, `Superscaffold_1`, `SUPER_1`
- `chr1`, `chromosome_1`, `Chr_1`
- `LG_1` (linkage groups)
- `scaffold_1_cov50` (coverage-annotated)
- `HiC_scaffold_1`, `Scaffold_1_RaGOO`
- `NC_000001.11`, `CM000001.1` (NCBI accessions)

This inconsistency makes automated analysis and cross-species comparisons difficult. Existing QC tools like QUAST report metrics but don't classify scaffolds. Scaffolding tools like LACHESIS create assemblies but don't help interpret existing ones.

## The Solution

ChromDetect uses multiple complementary strategies to identify chromosome-level scaffolds:

1. **Name-based detection** - Regex patterns for 15+ common naming conventions
2. **Size-based detection** - Large scaffolds are typically chromosomes
3. **N50-based detection** - Scaffolds contributing to N50 are typically chromosome-level
4. **Karyotype-informed detection** - Use known chromosome count to adjust classifications

## Installation

```bash
pip install chromdetect
```

Or install from source:

```bash
git clone https://github.com/verity-project/chromdetect.git
cd chromdetect
pip install -e .
```

## Quick Start

### Command Line

```bash
# Basic usage - get summary
chromdetect assembly.fasta

# Output JSON for programmatic use
chromdetect assembly.fasta --format json --output results.json

# Use karyotype information for better accuracy
chromdetect assembly.fasta --karyotype 24

# Export only chromosome-level scaffolds as TSV
chromdetect assembly.fasta --chromosomes-only --format tsv > chromosomes.tsv
```

### Python API

```python
from chromdetect import parse_fasta, classify_scaffolds

# Parse and classify
scaffolds = parse_fasta("assembly.fasta.gz")
results, stats = classify_scaffolds(scaffolds, expected_chromosomes=24)

# Print summary
print(f"Found {stats.chromosome_count} chromosomes")
print(f"Total assembly: {stats.total_length / 1e9:.2f} Gb")
print(f"N50: {stats.n50 / 1e6:.1f} Mb")

# Access individual scaffold classifications
for r in results:
    if r.classification == "chromosome":
        print(f"{r.name}: {r.length:,} bp (confidence: {r.confidence:.2f})")
```

## Output Formats

### Summary (default)

```
============================================================
CHROMDETECT ASSEMBLY ANALYSIS
============================================================

Total scaffolds:     1,234
Total length:        2,876,543,210 bp (2.88 Gb)
N50:                 45,678,901 bp (45.7 Mb)
N90:                 12,345,678 bp
Largest scaffold:    198,765,432 bp

Scaffold Classification:
  Chromosomes:       24 (2.85 Gb)
  Unlocalized:       15
  Unplaced:          1,195

Chromosome N50:      118,234,567 bp (118.2 Mb)
GC content:          41.2%
```

### JSON

```json
{
  "summary": {
    "total_scaffolds": 1234,
    "chromosome_count": 24,
    "n50": 45678901,
    ...
  },
  "scaffolds": [
    {
      "name": "chr1",
      "length": 198765432,
      "classification": "chromosome",
      "confidence": 0.95,
      "detection_method": "name_chr_explicit",
      "chromosome_id": "1"
    },
    ...
  ]
}
```

### TSV

```
name    length    classification    confidence    method    chromosome_id
chr1    198765432    chromosome    0.95    name_chr_explicit    1
chr2    175432198    chromosome    0.93    name_chr_explicit    2
...
```

## Options

| Option | Description |
|--------|-------------|
| `-f, --format` | Output format: `summary`, `json`, `tsv` (default: summary) |
| `-o, --output` | Write output to file instead of stdout |
| `-k, --karyotype` | Expected chromosome count for karyotype-informed detection |
| `-s, --min-size` | Minimum size (bp) to consider chromosome-level (default: 10Mb) |
| `-c, --chromosomes-only` | Only output chromosome-level scaffolds |
| `-q, --quiet` | Suppress progress messages |

## Supported Naming Conventions

ChromDetect recognizes these naming patterns (case-insensitive):

| Pattern | Examples | Method |
|---------|----------|--------|
| Explicit chromosome | `chr1`, `chromosome_X`, `Chr_MT` | `name_chr_explicit` |
| Super scaffold | `Super_scaffold_1`, `Superscaffold_X` | `name_super_scaffold` |
| SUPER | `SUPER_1`, `SUPER1` | `name_SUPER` |
| Linkage group | `LG1`, `LG_X` | `name_linkage_group` |
| NCBI RefSeq | `NC_000001.11` | `name_ncbi_refseq` |
| NCBI GenBank | `CM000001.1` | `name_ncbi_genbank` |
| HiC scaffold | `HiC_scaffold_1` | `name_hic_scaffold` |
| RaGOO | `Scaffold_1_RaGOO` | `name_ragoo` |
| Simple numeric | `1`, `X`, `MT` | `name_numeric` |

Patterns that indicate **unlocalized** scaffolds:
- `*_random`, `*_unloc*`, `chrUn_*`

Patterns that indicate **unplaced** scaffolds (contigs/fragments):
- `*_ctg*`, `*contig*`, `*_arrow_*`, `*_pilon*`, `*_hap*`

## How It Works

ChromDetect combines name-based and size-based detection with these priority rules:

1. **Strong name match (confidence ≥ 0.8)** takes priority
2. **Large scaffold + weak name match** = chromosome with boosted confidence
3. **Large scaffold + no name match** = chromosome with reduced confidence
4. **Small scaffold** = unplaced regardless of name

When `--karyotype` is provided:
- If too many candidates: demote lowest-confidence chromosomes
- If too few candidates: promote largest unplaced scaffolds

## Use Cases

### VGP Assembly Validation

```bash
# Validate a VGP curated assembly
chromdetect species.pri.cur.fasta.gz --karyotype 24 --format json
```

### Cross-Species Comparison

```python
from chromdetect import parse_fasta, classify_scaffolds

species_files = ["human.fa", "mouse.fa", "zebrafish.fa"]
karyotypes = [23, 20, 25]

for fasta, n_chr in zip(species_files, karyotypes):
    scaffolds = parse_fasta(fasta)
    results, stats = classify_scaffolds(scaffolds, expected_chromosomes=n_chr)
    print(f"{fasta}: {stats.chromosome_count} chromosomes detected")
```

### Pipeline Integration

```bash
# As part of assembly QC pipeline
chromdetect assembly.fasta --format json | jq '.summary.chromosome_count'
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Adding New Patterns

To add support for a new naming convention:

1. Add the regex pattern to `chromdetect/patterns.py`
2. Include a descriptive method name
3. Ensure the pattern captures chromosome ID in group 1
4. Add tests in `tests/test_patterns.py`

Example:

```python
# In patterns.py
CHROMOSOME_PATTERNS.append(
    (r'^MyConvention_(\d+)$', 'my_convention'),
)
```

## Citation

If you use ChromDetect in your research, please cite:

```
ChromDetect: Chromosome-level scaffold detection for genome assemblies
https://github.com/verity-project/chromdetect
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Related Projects

- [QUAST](https://github.com/ablab/quast) - Quality assessment tool for genome assemblies
- [Verity](https://github.com/verity-project/verity) - Hi-C-based assembly validation framework
