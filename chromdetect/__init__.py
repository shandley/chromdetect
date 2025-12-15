"""
ChromDetect - Chromosome-Level Scaffold Detection Tool

A tool to identify chromosome-level scaffolds in genome assemblies
using multiple detection strategies to handle inconsistent naming conventions.
"""

from chromdetect.assembly_report import (
    AssemblyReport,
    AssemblyReportEntry,
    apply_assembly_report,
    parse_assembly_report,
)
from chromdetect.compare import (
    ComparisonResult,
    compare_assemblies,
    compare_fasta_files,
    format_comparison_summary,
    format_comparison_tsv,
)
from chromdetect.core import (
    AssemblyStats,
    ScaffoldInfo,
    classify_fasta,
    classify_scaffolds,
    detect_by_name,
    detect_by_size,
    format_bed,
    format_gff,
    parse_fasta,
    parse_fasta_from_handle,
    write_fasta,
)
from chromdetect.html_report import generate_html_report
from chromdetect.patterns import (
    CHROMOSOME_PATTERNS,
    FRAGMENT_PATTERNS,
    UNLOCALIZED_PATTERNS,
    load_custom_patterns,
    merge_patterns,
)

__version__ = "0.5.0"
__all__ = [
    # Core classification
    "ScaffoldInfo",
    "AssemblyStats",
    "classify_scaffolds",
    "classify_fasta",  # Convenience function
    "parse_fasta",
    "parse_fasta_from_handle",
    "detect_by_name",
    "detect_by_size",
    "write_fasta",
    "format_bed",
    "format_gff",
    # Patterns
    "CHROMOSOME_PATTERNS",
    "UNLOCALIZED_PATTERNS",
    "FRAGMENT_PATTERNS",
    "load_custom_patterns",
    "merge_patterns",
    # Assembly reports
    "AssemblyReport",
    "AssemblyReportEntry",
    "parse_assembly_report",
    "apply_assembly_report",
    # Comparison
    "ComparisonResult",
    "compare_assemblies",
    "compare_fasta_files",  # Convenience function
    "format_comparison_summary",
    "format_comparison_tsv",
    # HTML report
    "generate_html_report",
]
