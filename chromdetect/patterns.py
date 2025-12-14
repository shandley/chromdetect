"""
Chromosome naming pattern definitions for ChromDetect.

This module contains regex patterns for detecting chromosome-level scaffolds
across different naming conventions used by various genome assembly pipelines.

Adding New Patterns:
    To add support for a new naming convention:
    1. Add the regex pattern to CHROMOSOME_PATTERNS
    2. Include a descriptive method name
    3. Ensure the pattern captures the chromosome ID in group 1
    4. Add tests in tests/test_patterns.py
"""

import re
from typing import Pattern

# Common chromosome naming patterns (case-insensitive by default)
# Format: (pattern, method_name)
# The pattern should capture chromosome ID (number, X, Y, etc.) in group 1 when possible
CHROMOSOME_PATTERNS: list[tuple[str, str]] = [
    # Explicit chromosome names (most reliable)
    (r'^chr(?:omosome)?[_\-\s]?(\d+|[XYZWAB]|MT?|Un)$', 'chr_explicit'),

    # VGP-style naming conventions
    (r'^super[_\-\s]?scaffold[_\-\s]?(\d+|[XYZWAB])$', 'super_scaffold'),
    (r'^superscaffold[_\-\s]?(\d+|[XYZWAB])$', 'superscaffold'),
    (r'^SUPER[_\-\s]?(\d+|[XYZWAB])$', 'SUPER'),

    # Linkage group naming (common in fish, plants)
    (r'^LG[_\-\s]?(\d+|[XYZWAB])$', 'linkage_group'),

    # GenBank/RefSeq accession patterns (assembled chromosomes)
    (r'^NC_\d+\.\d+$', 'ncbi_refseq'),
    (r'^CM\d+\.\d+$', 'ncbi_genbank'),

    # Simple numeric chromosome names
    (r'^(?:chr)?(\d+|[XYZWAB])$', 'numeric'),

    # HiC scaffolder output patterns
    (r'^HiC_scaffold_(\d+)$', 'hic_scaffold'),
    (r'^Scaffold_(\d+)_RaGOO$', 'ragoo'),

    # Plant genome conventions
    (r'^Pt(\d+)$', 'plant_chromosome'),
    (r'^Gm(\d+)$', 'soybean_chromosome'),

    # Additional VGP variations
    (r'^Super-Scaffold_(\d+)$', 'super_scaffold_hyphen'),
    (r'^scaffold_(\d+)_cov\d+$', 'scaffold_cov'),
]

# Patterns suggesting unlocalized/random scaffolds (chromosome-associated but not placed)
UNLOCALIZED_PATTERNS: list[str] = [
    r'random',
    r'unloc',
    r'unplaced',
    r'_un_',
    r'chrUn',
    r'scaffold.*unloc',
    r'_hap\d+_unloc',
]

# Patterns suggesting contigs/fragments (not chromosome-level)
FRAGMENT_PATTERNS: list[str] = [
    r'ctg\d*$',
    r'contig',
    r'_arrow_',
    r'_pilon',
    r'fragment',
    r'_hap\d',
    r'_alt$',
    r'_patch$',
    r'debris',
]

# Sex chromosome identifiers
SEX_CHROMOSOMES: set[str] = {'X', 'Y', 'Z', 'W', 'x', 'y', 'z', 'w'}

# Mitochondrial identifiers
MITOCHONDRIAL_IDS: set[str] = {'M', 'MT', 'Mt', 'mt', 'mito', 'mitochondrion'}


def compile_patterns() -> list[tuple[Pattern[str], str]]:
    """Compile chromosome patterns for efficient matching."""
    return [(re.compile(p, re.IGNORECASE), name) for p, name in CHROMOSOME_PATTERNS]


def compile_exclusion_patterns() -> tuple[list[Pattern[str]], list[Pattern[str]]]:
    """Compile exclusion patterns (unlocalized, fragments)."""
    unloc = [re.compile(p, re.IGNORECASE) for p in UNLOCALIZED_PATTERNS]
    frag = [re.compile(p, re.IGNORECASE) for p in FRAGMENT_PATTERNS]
    return unloc, frag


# Pre-compiled patterns for performance
COMPILED_CHROMOSOME_PATTERNS = compile_patterns()
COMPILED_UNLOCALIZED, COMPILED_FRAGMENT = compile_exclusion_patterns()
