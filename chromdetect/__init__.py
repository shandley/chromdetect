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
    format_bed,
    format_gff,
    parse_fasta,
    parse_fasta_from_handle,
    write_fasta,
)
from chromdetect.patterns import (
    CHROMOSOME_PATTERNS,
    FRAGMENT_PATTERNS,
    UNLOCALIZED_PATTERNS,
)

__version__ = "0.3.0"
__all__ = [
    "ScaffoldInfo",
    "AssemblyStats",
    "classify_scaffolds",
    "parse_fasta",
    "parse_fasta_from_handle",
    "detect_by_name",
    "detect_by_size",
    "write_fasta",
    "format_bed",
    "format_gff",
    "CHROMOSOME_PATTERNS",
    "UNLOCALIZED_PATTERNS",
    "FRAGMENT_PATTERNS",
]
