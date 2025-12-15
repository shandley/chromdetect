# ChromDetect Example Data

This directory contains synthetic genome assemblies for testing and demonstration purposes.

## Included Files

### `synthetic_assembly.fasta`

A small (~5 KB) synthetic assembly demonstrating various scaffold naming conventions:

| Scaffold | Type | Naming Convention |
|----------|------|-------------------|
| chr1 | Chromosome | Explicit chromosome |
| chr2 | Chromosome | Explicit chromosome |
| chrX | Sex chromosome | Explicit chromosome |
| chrMT | Mitochondrial | Explicit chromosome |
| Super_scaffold_1 | Chromosome | VGP-style |
| HiC_scaffold_1 | Chromosome | Hi-C scaffolder |
| LG1 | Chromosome | Linkage group |
| NC_000001.11 | Chromosome | NCBI RefSeq |
| scaffold_1_RaGOO | Chromosome | RaGOO scaffolder |
| chr1_random | Unlocalized | Random scaffold |
| chrUn_scaffold1 | Unplaced | Unknown chromosome |
| contig_001 | Unplaced | Contig |
| scaffold_arrow_1 | Unplaced | Arrow-polished |

### `synthetic_assembly_v2.fasta`

An "improved" version of the assembly for testing comparison mode:
- chr1: Extended
- chr2: Extended
- chrX: Extended
- chr3: New chromosome added
- Some small scaffolds removed (merged)

## Quick Start

```bash
# Basic analysis
chromdetect examples/synthetic_assembly.fasta

# Compare the two versions
chromdetect examples/synthetic_assembly.fasta --compare examples/synthetic_assembly_v2.fasta

# Generate HTML report
chromdetect examples/synthetic_assembly.fasta --format html -o report.html

# Expected karyotype (10 chromosomes)
chromdetect examples/synthetic_assembly.fasta --karyotype 10
```

## Expected Output

Running `chromdetect examples/synthetic_assembly.fasta` should show:
- **Chromosomes detected:** 9 (chr1, chr2, chrX, chrMT, Super_scaffold_1, HiC_scaffold_1, LG1, NC_000001.11, scaffold_1_RaGOO)
- **Unlocalized:** 2 (chr1_random, chrUn_scaffold1)
- **Unplaced:** 2 (contig_001, scaffold_arrow_1)

## Real-World Test Data

For testing with real genome assemblies, see the main README for download instructions for:
- *Saccharomyces cerevisiae* S288C (~12 Mb) - Yeast reference genome
- *Caenorhabditis elegans* (~100 Mb) - Nematode reference genome
- *Arabidopsis thaliana* (~135 Mb) - Plant reference genome
