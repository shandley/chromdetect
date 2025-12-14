"""
ChromDetect - Chromosome-Level Scaffold Detection Tool

A tool to identify chromosome-level scaffolds in genome assemblies
using multiple detection strategies to handle inconsistent naming conventions.
"""

from chromdetect.core import (
    AssemblyStats,
    ScaffoldInfo,
    classify_scaffolds,
    detect_by_name,
    detect_by_size,
    parse_fasta,
    parse_fasta_from_handle,
)
from chromdetect.patterns import (
    CHROMOSOME_PATTERNS,
    FRAGMENT_PATTERNS,
    UNLOCALIZED_PATTERNS,
)

__version__ = "0.2.0"
__all__ = [
    "ScaffoldInfo",
    "AssemblyStats",
    "classify_scaffolds",
    "parse_fasta",
    "parse_fasta_from_handle",
    "detect_by_name",
    "detect_by_size",
    "CHROMOSOME_PATTERNS",
    "UNLOCALIZED_PATTERNS",
    "FRAGMENT_PATTERNS",
]
