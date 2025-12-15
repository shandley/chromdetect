"""
Assembly comparison module for ChromDetect.

This module provides functionality to compare two genome assemblies
side-by-side, analyzing differences in scaffold classification,
statistics, and completeness.
"""

from __future__ import annotations

from dataclasses import dataclass

from chromdetect.core import AssemblyStats, ScaffoldInfo


@dataclass
class ComparisonResult:
    """Result of comparing two assemblies.

    Attributes:
        assembly1_name: Name of first assembly
        assembly2_name: Name of second assembly
        stats1: Statistics for first assembly
        stats2: Statistics for second assembly
        shared_chromosomes: Chromosomes found in both assemblies
        unique_to_1: Chromosomes only in first assembly
        unique_to_2: Chromosomes only in second assembly
        size_differences: Dict mapping chromosome ID to size difference
        classification_changes: Scaffolds with different classifications
    """

    assembly1_name: str
    assembly2_name: str
    stats1: AssemblyStats
    stats2: AssemblyStats
    shared_chromosomes: list[str]
    unique_to_1: list[str]
    unique_to_2: list[str]
    size_differences: dict[str, int]
    classification_changes: list[tuple[str, str, str]]  # (name, class1, class2)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "assembly1_name": self.assembly1_name,
            "assembly2_name": self.assembly2_name,
            "stats1": self.stats1.to_dict(),
            "stats2": self.stats2.to_dict(),
            "shared_chromosomes": self.shared_chromosomes,
            "unique_to_1": self.unique_to_1,
            "unique_to_2": self.unique_to_2,
            "size_differences": self.size_differences,
            "classification_changes": [
                {"name": name, "classification1": c1, "classification2": c2}
                for name, c1, c2 in self.classification_changes
            ],
            "summary": self.summary(),
        }

    def summary(self) -> dict:
        """Generate summary statistics for the comparison."""
        return {
            "total_shared_chromosomes": len(self.shared_chromosomes),
            "unique_to_1": len(self.unique_to_1),
            "unique_to_2": len(self.unique_to_2),
            "classification_changes_count": len(self.classification_changes),
            "n50_difference": self.stats2.n50 - self.stats1.n50,
            "chromosome_count_difference": (
                self.stats2.chromosome_count - self.stats1.chromosome_count
            ),
            "total_length_difference": (
                self.stats2.total_length - self.stats1.total_length
            ),
        }


def compare_assemblies(
    results1: list[ScaffoldInfo],
    stats1: AssemblyStats,
    results2: list[ScaffoldInfo],
    stats2: AssemblyStats,
    assembly1_name: str = "Assembly 1",
    assembly2_name: str = "Assembly 2",
) -> ComparisonResult:
    """
    Compare two genome assemblies.

    Args:
        results1: ScaffoldInfo list from first assembly
        stats1: AssemblyStats from first assembly
        results2: ScaffoldInfo list from second assembly
        stats2: AssemblyStats from second assembly
        assembly1_name: Display name for first assembly
        assembly2_name: Display name for second assembly

    Returns:
        ComparisonResult with detailed comparison information
    """
    # Get chromosome scaffolds from each assembly
    chr1 = {r.name: r for r in results1 if r.classification == "chromosome"}
    chr2 = {r.name: r for r in results2 if r.classification == "chromosome"}

    # Find shared and unique chromosomes
    shared_names = set(chr1.keys()) & set(chr2.keys())
    unique_to_1 = sorted(set(chr1.keys()) - set(chr2.keys()))
    unique_to_2 = sorted(set(chr2.keys()) - set(chr1.keys()))
    shared_chromosomes = sorted(shared_names)

    # Calculate size differences for shared chromosomes
    size_differences = {}
    for name in shared_chromosomes:
        size_diff = chr2[name].length - chr1[name].length
        if size_diff != 0:
            size_differences[name] = size_diff

    # Find scaffolds with different classifications
    # Build lookup by scaffold name across all scaffolds
    all_names_1 = {r.name: r for r in results1}
    all_names_2 = {r.name: r for r in results2}

    classification_changes = []
    # Check scaffolds present in both
    common_scaffolds = set(all_names_1.keys()) & set(all_names_2.keys())
    for name in sorted(common_scaffolds):
        r1 = all_names_1[name]
        r2 = all_names_2[name]
        if r1.classification != r2.classification:
            classification_changes.append((name, r1.classification, r2.classification))

    return ComparisonResult(
        assembly1_name=assembly1_name,
        assembly2_name=assembly2_name,
        stats1=stats1,
        stats2=stats2,
        shared_chromosomes=shared_chromosomes,
        unique_to_1=unique_to_1,
        unique_to_2=unique_to_2,
        size_differences=size_differences,
        classification_changes=classification_changes,
    )


def format_comparison_summary(comparison: ComparisonResult) -> str:
    """
    Format comparison result as human-readable summary.

    Args:
        comparison: ComparisonResult from compare_assemblies

    Returns:
        Formatted string for display
    """
    lines = [
        "=" * 70,
        "CHROMDETECT ASSEMBLY COMPARISON",
        "=" * 70,
        "",
        f"Assembly 1: {comparison.assembly1_name}",
        f"Assembly 2: {comparison.assembly2_name}",
        "",
        "-" * 70,
        "STATISTICS COMPARISON",
        "-" * 70,
        "",
        f"{'Metric':<30} {'Assembly 1':>18} {'Assembly 2':>18}",
        "-" * 70,
    ]

    s1 = comparison.stats1
    s2 = comparison.stats2

    # Format statistics comparison
    lines.append(
        f"{'Total scaffolds':<30} {s1.total_scaffolds:>18,} {s2.total_scaffolds:>18,}"
    )
    lines.append(
        f"{'Total length (bp)':<30} {s1.total_length:>18,} {s2.total_length:>18,}"
    )
    lines.append(f"{'N50 (bp)':<30} {s1.n50:>18,} {s2.n50:>18,}")
    lines.append(f"{'N90 (bp)':<30} {s1.n90:>18,} {s2.n90:>18,}")
    lines.append(
        f"{'Chromosome count':<30} {s1.chromosome_count:>18} {s2.chromosome_count:>18}"
    )
    lines.append(
        f"{'Chromosome length (bp)':<30} {s1.chromosome_length:>18,} {s2.chromosome_length:>18,}"
    )
    lines.append(
        f"{'Chromosome N50 (bp)':<30} {s1.chromosome_n50:>18,} {s2.chromosome_n50:>18,}"
    )
    lines.append(
        f"{'Unlocalized count':<30} {s1.unlocalized_count:>18} {s2.unlocalized_count:>18}"
    )
    lines.append(
        f"{'Unplaced count':<30} {s1.unplaced_count:>18} {s2.unplaced_count:>18}"
    )

    # GC content
    gc1 = f"{s1.gc_content * 100:.2f}%" if s1.gc_content else "N/A"
    gc2 = f"{s2.gc_content * 100:.2f}%" if s2.gc_content else "N/A"
    lines.append(f"{'GC content':<30} {gc1:>18} {gc2:>18}")

    # Chromosome comparison section
    lines.extend(
        [
            "",
            "-" * 70,
            "CHROMOSOME COMPARISON",
            "-" * 70,
            "",
            f"Shared chromosomes:      {len(comparison.shared_chromosomes)}",
            f"Unique to Assembly 1:    {len(comparison.unique_to_1)}",
            f"Unique to Assembly 2:    {len(comparison.unique_to_2)}",
        ]
    )

    # Show unique chromosomes
    if comparison.unique_to_1:
        lines.append("")
        lines.append(f"Chromosomes only in {comparison.assembly1_name}:")
        for name in comparison.unique_to_1[:10]:  # Show first 10
            lines.append(f"  - {name}")
        if len(comparison.unique_to_1) > 10:
            lines.append(f"  ... and {len(comparison.unique_to_1) - 10} more")

    if comparison.unique_to_2:
        lines.append("")
        lines.append(f"Chromosomes only in {comparison.assembly2_name}:")
        for name in comparison.unique_to_2[:10]:  # Show first 10
            lines.append(f"  - {name}")
        if len(comparison.unique_to_2) > 10:
            lines.append(f"  ... and {len(comparison.unique_to_2) - 10} more")

    # Size differences
    if comparison.size_differences:
        lines.extend(
            [
                "",
                "-" * 70,
                "SIZE DIFFERENCES (shared chromosomes)",
                "-" * 70,
                "",
            ]
        )

        # Sort by absolute difference
        sorted_diffs = sorted(
            comparison.size_differences.items(), key=lambda x: abs(x[1]), reverse=True
        )
        for name, diff in sorted_diffs[:10]:
            sign = "+" if diff > 0 else ""
            lines.append(f"  {name:<30} {sign}{diff:>12,} bp")
        if len(sorted_diffs) > 10:
            lines.append(f"  ... and {len(sorted_diffs) - 10} more")

    # Classification changes
    if comparison.classification_changes:
        lines.extend(
            [
                "",
                "-" * 70,
                "CLASSIFICATION CHANGES",
                "-" * 70,
                "",
            ]
        )
        for name, c1, c2 in comparison.classification_changes[:10]:
            lines.append(f"  {name:<30} {c1:<15} -> {c2:<15}")
        if len(comparison.classification_changes) > 10:
            lines.append(f"  ... and {len(comparison.classification_changes) - 10} more")

    # Summary
    summary = comparison.summary()
    lines.extend(
        [
            "",
            "-" * 70,
            "SUMMARY",
            "-" * 70,
            "",
        ]
    )

    # N50 improvement/degradation
    n50_diff = summary["n50_difference"]
    if n50_diff > 0:
        lines.append(f"N50 improved by {n50_diff:,} bp in {comparison.assembly2_name}")
    elif n50_diff < 0:
        lines.append(
            f"N50 decreased by {abs(n50_diff):,} bp in {comparison.assembly2_name}"
        )
    else:
        lines.append("N50 unchanged between assemblies")

    lines.append("")

    return "\n".join(lines)


def format_comparison_tsv(comparison: ComparisonResult) -> str:
    """
    Format comparison result as TSV.

    Args:
        comparison: ComparisonResult from compare_assemblies

    Returns:
        TSV-formatted string
    """
    lines = [
        "metric\tassembly1\tassembly2\tdifference",
    ]

    s1 = comparison.stats1
    s2 = comparison.stats2

    lines.append(
        f"total_scaffolds\t{s1.total_scaffolds}\t{s2.total_scaffolds}\t"
        f"{s2.total_scaffolds - s1.total_scaffolds}"
    )
    lines.append(
        f"total_length\t{s1.total_length}\t{s2.total_length}\t"
        f"{s2.total_length - s1.total_length}"
    )
    lines.append(f"n50\t{s1.n50}\t{s2.n50}\t{s2.n50 - s1.n50}")
    lines.append(f"n90\t{s1.n90}\t{s2.n90}\t{s2.n90 - s1.n90}")
    lines.append(
        f"chromosome_count\t{s1.chromosome_count}\t{s2.chromosome_count}\t"
        f"{s2.chromosome_count - s1.chromosome_count}"
    )
    lines.append(
        f"chromosome_length\t{s1.chromosome_length}\t{s2.chromosome_length}\t"
        f"{s2.chromosome_length - s1.chromosome_length}"
    )
    lines.append(
        f"chromosome_n50\t{s1.chromosome_n50}\t{s2.chromosome_n50}\t"
        f"{s2.chromosome_n50 - s1.chromosome_n50}"
    )

    gc1 = s1.gc_content or 0
    gc2 = s2.gc_content or 0
    lines.append(f"gc_content\t{gc1:.4f}\t{gc2:.4f}\t{gc2 - gc1:.4f}")

    lines.append(
        f"shared_chromosomes\t{len(comparison.shared_chromosomes)}\t"
        f"{len(comparison.shared_chromosomes)}\t0"
    )
    lines.append(
        f"unique_chromosomes\t{len(comparison.unique_to_1)}\t"
        f"{len(comparison.unique_to_2)}\t"
        f"{len(comparison.unique_to_2) - len(comparison.unique_to_1)}"
    )

    return "\n".join(lines)


def compare_fasta_files(
    fasta1_path: str,
    fasta2_path: str,
) -> ComparisonResult:
    """
    Convenience function to compare two FASTA files directly.

    This is a high-level API that handles parsing and classification
    internally, making it easy to compare two assemblies with a single call.

    Args:
        fasta1_path: Path to first FASTA file (can be gzipped)
        fasta2_path: Path to second FASTA file (can be gzipped)

    Returns:
        ComparisonResult with detailed comparison information

    Example:
        >>> result = compare_fasta_files("assembly_v1.fasta", "assembly_v2.fasta")
        >>> print(f"N50 improved by {result.summary()['n50_difference']} bp")
    """
    from pathlib import Path

    from chromdetect.core import classify_scaffolds, parse_fasta

    # Parse and classify first assembly
    scaffolds1 = parse_fasta(fasta1_path)
    results1, stats1 = classify_scaffolds(scaffolds1)

    # Parse and classify second assembly
    scaffolds2 = parse_fasta(fasta2_path)
    results2, stats2 = classify_scaffolds(scaffolds2)

    # Extract assembly names from paths
    name1 = Path(fasta1_path).stem.replace(".fasta", "").replace(".fa", "")
    name2 = Path(fasta2_path).stem.replace(".fasta", "").replace(".fa", "")

    return compare_assemblies(
        results1, stats1,
        results2, stats2,
        assembly1_name=name1,
        assembly2_name=name2,
    )
