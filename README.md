<p align="center">
  <img src="chromdetect_hex_logo.jpeg" alt="ChromDetect Logo" width="200"/>
</p>

<h1 align="center">ChromDetect</h1>

<p align="center">
  <a href="https://pypi.org/project/chromdetect/"><img src="https://img.shields.io/pypi/v/chromdetect.svg" alt="PyPI version"></a>
  <a href="https://pypi.org/project/chromdetect/"><img src="https://img.shields.io/pypi/pyversions/chromdetect.svg" alt="Python versions"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://github.com/shandley/chromdetect/actions/workflows/test.yml"><img src="https://github.com/shandley/chromdetect/actions/workflows/test.yml/badge.svg" alt="Tests"></a>
</p>

<p align="center">
  <strong>Detect chromosome-level scaffolds in genome assemblies with inconsistent naming conventions.</strong>
</p>

---

## The Problem

Genome assemblies use wildly inconsistent naming conventions for chromosome-level scaffolds:

- `Super_scaffold_1`, `Superscaffold_1`, `SUPER_1`
- `chr1`, `chromosome_1`, `Chr_1`
- `LG_1` (linkage groups)
- `scaffold_1_cov50` (coverage-annotated)
- `HiC_scaffold_1`, `Scaffold_1_RaGOO`
- `NC_000001.11`, `CM000001.1` (NCBI accessions)

This inconsistency makes automated analysis and cross-species comparisons difficult. Existing QC tools like QUAST report metrics but don't classify scaffolds. Scaffolding tools like LACHESIS create assemblies but don't help interpret existing ones.

## Why ChromDetect?

| Feature | QUAST | assembly-stats | gfastats | **ChromDetect** |
|---------|-------|----------------|----------|-----------------|
| N50/N90 statistics | ✅ | ✅ | ✅ | ✅ |
| Scaffold classification | ❌ | ❌ | ❌ | ✅ |
| Pattern-based detection | ❌ | ❌ | ❌ | ✅ |
| Size-based detection | ❌ | ❌ | ❌ | ✅ |
| Karyotype-aware | ❌ | ❌ | ❌ | ✅ |
| Telomere detection (T2T) | ❌ | ❌ | ❌ | ✅ |
| Quality scoring | ❌ | ❌ | ❌ | ✅ |
| NCBI report integration | ❌ | ❌ | ❌ | ✅ |
| Multiple output formats | ✅ | ❌ | ✅ | ✅ |
| Zero dependencies | ❌ | ✅ | ❌ | ✅ |

ChromDetect fills a gap in the genomics toolkit: **automatically identifying which scaffolds represent chromosomes** rather than just reporting assembly statistics.

## The Solution

ChromDetect uses multiple complementary strategies to identify chromosome-level scaffolds:

1. **Name-based detection** - Regex patterns for 15+ common naming conventions
2. **Size-based detection** - Large scaffolds are typically chromosomes
3. **N50-based detection** - Scaffolds contributing to N50 are typically chromosome-level
4. **Karyotype-informed detection** - Use known chromosome count to adjust classifications
5. **Telomere detection** - Identify T2T chromosomes by detecting telomeric repeats at scaffold ends
6. **NCBI assembly report** - Use official NCBI reports for authoritative classification

## Installation

```bash
pip install chromdetect
```

Or install from source:

```bash
git clone https://github.com/shandley/chromdetect.git
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

# Export as BED or GFF format for pipeline integration
chromdetect assembly.fasta --format bed > scaffolds.bed
chromdetect assembly.fasta --format gff > scaffolds.gff

# Extract chromosome sequences to a new FASTA file
chromdetect assembly.fasta --extract-chromosomes chromosomes.fasta

# Batch process multiple assemblies
chromdetect --batch assemblies_dir/ --output results_dir/

# Detect telomeres at scaffold ends (T2T detection)
chromdetect assembly.fasta --detect-telomeres

# Use custom naming patterns
chromdetect assembly.fasta --patterns custom_patterns.yaml

# Use NCBI assembly report for accurate classification
chromdetect assembly.fasta --assembly-report GCF_000001405.assembly_report.txt
```

### Python API

```python
from chromdetect import (
    parse_fasta, classify_scaffolds, write_fasta, format_bed, format_gff,
    detect_telomere, parse_assembly_report, calculate_quality_score
)

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

# Export to BED or GFF format
bed_output = format_bed(results)
gff_output = format_gff(results)

# Extract chromosome sequences (requires full sequence parsing)
scaffolds = parse_fasta("assembly.fasta", keep_full_sequence=True)
results, stats = classify_scaffolds(scaffolds)
chr_seqs = [(name, seq) for name, length, seq in scaffolds
            if any(r.name == name and r.classification == "chromosome" for r in results)]
write_fasta(chr_seqs, "chromosomes.fasta")

# Enable telomere detection for T2T identification
scaffolds = parse_fasta("assembly.fasta", keep_full_sequence=True)
results, stats = classify_scaffolds(scaffolds, detect_telomeres=True)
for r in results:
    if r.classification == "chromosome":
        status = "T2T" if r.is_t2t else ("Telomere" if r.has_telomere else "")
        print(f"{r.name}: {r.length:,} bp {status}")

# Use NCBI assembly report for authoritative classification
report = parse_assembly_report("assembly_report.txt")
results, stats = classify_scaffolds(scaffolds, assembly_report=report)
print(f"Quality score: {stats.quality_score:.3f}")
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

Telomere Detection:
  With telomeres:    22
  T2T chromosomes:   18

Quality Score:       0.847
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

### BED

Standard BED6 format for integration with bedtools, IGV, and other genomics tools:

```
chr1    0    198765432    chromosome    950    .
chr2    0    175432198    chromosome    930    .
...
```

### GFF3

GFF3 format with classification metadata in attributes:

```
##gff-version 3
chr1    chromdetect    chromosome    1    198765432    0.950    .    .    ID=chr1;Name=chr1;classification=chromosome;detection_method=name_chr_explicit;chromosome_id=1
...
```

## Options

| Option | Description |
|--------|-------------|
| `-f, --format` | Output format: `summary`, `json`, `tsv`, `bed`, `gff` (default: summary) |
| `-o, --output` | Write output to file instead of stdout |
| `-k, --karyotype` | Expected chromosome count for karyotype-informed detection |
| `-s, --min-size` | Minimum size (bp) to consider chromosome-level (default: 10Mb) |
| `-c, --chromosomes-only` | Only output chromosome-level scaffolds |
| `--extract-chromosomes` | Extract chromosome sequences to a FASTA file |
| `--batch` | Process all FASTA files in a directory |
| `--detect-telomeres` | Detect telomeric repeats at scaffold ends for T2T identification |
| `--patterns` | Custom patterns file (YAML or JSON) for scaffold name matching |
| `--assembly-report` | NCBI assembly report file for authoritative classification |
| `--min-confidence` | Minimum confidence threshold (0.0-1.0) to include scaffolds |
| `--min-length` | Minimum scaffold length (bp) to include in output |
| `-q, --quiet` | Suppress progress messages |
| `-v, --verbose` | Show detailed processing information |

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

### Telomere Detection

When `--detect-telomeres` is enabled, ChromDetect searches scaffold ends for telomeric repeat motifs:

| Organism Type | 5' Motif | 3' Motif |
|--------------|----------|----------|
| Vertebrates | CCCTAA | TTAGGG |
| Plants (Arabidopsis) | CCCTAAA | TTTAGGG |
| Insects (Bombyx) | CCTAA | TTAGG |
| Nematodes | GCCTAA | TTAGGC |

Scaffolds with telomeres at both ends are marked as **T2T** (telomere-to-telomere), indicating complete chromosome assembly.

### Quality Score

ChromDetect calculates an overall assembly quality score (0.0-1.0) based on:

- **Classification confidence** (30%) - Average confidence of chromosome classifications
- **Karyotype completeness** (25%) - How close to expected chromosome count
- **Telomere completeness** (25%) - Proportion of T2T chromosomes (when `--detect-telomeres` used)
- **Size consistency** (20%) - Chromosome length relative to total assembly

A score > 0.8 indicates a high-quality chromosome-level assembly.

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

# Export scaffold regions in BED format for downstream analysis
chromdetect assembly.fasta --format bed --chromosomes-only > chromosomes.bed
bedtools getfasta -fi assembly.fasta -bed chromosomes.bed -fo chr_regions.fa
```

### Batch Processing

```bash
# Process all assemblies in a directory
chromdetect --batch assemblies/ --format json --output results/

# This creates:
# - results/assembly1.json
# - results/assembly2.json
# - ...
# - results/batch_summary.tsv  (overview of all assemblies)
```

### Extract Chromosome Sequences

```bash
# Extract only chromosome-level sequences to a new FASTA
chromdetect assembly.fasta --extract-chromosomes chromosomes.fasta

# Combine with other options
chromdetect assembly.fasta \
    --karyotype 24 \
    --extract-chromosomes chromosomes.fasta \
    --format json --output report.json
```

### T2T Assembly Validation

```bash
# Detect telomeres at scaffold ends to identify T2T chromosomes
chromdetect assembly.fasta --detect-telomeres --format summary

# Full T2T analysis with quality scoring
chromdetect assembly.fasta \
    --detect-telomeres \
    --karyotype 23 \
    --format json --output t2t_report.json
```

### Using NCBI Assembly Reports

```bash
# Download an assembly report from NCBI
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40

# Use it for authoritative scaffold classification
chromdetect GRCh38.fasta --assembly-report GCF_000001405.40_GRCh38.p14_assembly_report.txt
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

### Using Custom Patterns

You can also use custom patterns without modifying the source code:

```yaml
# custom_patterns.yaml
chromosome_patterns:
  - pattern: "^MyScaffold_(\\d+)$"
    name: "my_scaffold"
  - pattern: "^CustomChr_(\\d+)$"
    name: "custom_chr"
unlocalized_patterns:
  - my_random
fragment_patterns:
  - my_contig
```

```bash
chromdetect assembly.fasta --patterns custom_patterns.yaml
```

## Citation

If you use ChromDetect in your research, please cite:

```
ChromDetect: Chromosome-level scaffold detection for genome assemblies
https://github.com/shandley/chromdetect
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Related Projects

- [QUAST](https://github.com/ablab/quast) - Quality assessment tool for genome assemblies
- [Verity](https://github.com/shandley/verity) - Hi-C-based assembly validation framework
