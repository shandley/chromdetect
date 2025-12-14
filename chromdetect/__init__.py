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
from chromdetect.core import (
    AssemblyStats,
    ScaffoldInfo,
    calculate_quality_score,
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
    load_custom_patterns,
    merge_patterns,
)
from chromdetect.telomere import (
    TELOMERE_MOTIFS,
    TelomereResult,
    detect_telomere,
    detect_telomeres_batch,
    get_telomere_summary,
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
    "calculate_quality_score",
    "detect_telomere",
    "detect_telomeres_batch",
    "get_telomere_summary",
    "TelomereResult",
    "TELOMERE_MOTIFS",
    "CHROMOSOME_PATTERNS",
    "UNLOCALIZED_PATTERNS",
    "FRAGMENT_PATTERNS",
    "load_custom_patterns",
    "merge_patterns",
    "AssemblyReport",
    "AssemblyReportEntry",
    "parse_assembly_report",
    "apply_assembly_report",
]
