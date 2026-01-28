"""
Assembly version tracking for ChromDetect.

This module provides functionality to compare two versions of a genome assembly,
tracking scaffold changes, promotions/demotions, and overall metric progression.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from chromdetect.core import (
    AssemblyStats,
    ScaffoldInfo,
    classify_scaffolds,
    parse_fasta,
)


class ChangeType(Enum):
    """Type of change between assembly versions."""

    UNCHANGED = "unchanged"  # Same name, same classification
    PROMOTED = "promoted"  # Moved to higher classification (e.g., unplaced → chromosome)
    DEMOTED = "demoted"  # Moved to lower classification (e.g., chromosome → unplaced)
    RESIZED = "resized"  # Significant size change (potential merge/split)
    RENAMED = "renamed"  # Name changed but likely same scaffold
    ADDED = "added"  # New scaffold in newer version
    REMOVED = "removed"  # Scaffold missing in newer version
    RECLASSIFIED = "reclassified"  # Same classification level, different type


# Classification hierarchy for promotion/demotion detection
CLASSIFICATION_RANK = {
    "chromosome": 3,
    "unlocalized": 2,
    "unplaced": 1,
    "other": 0,
}


@dataclass
class ScaffoldChange:
    """Represents a change to a single scaffold between versions.

    Attributes:
        name_v1: Scaffold name in version 1 (None if added)
        name_v2: Scaffold name in version 2 (None if removed)
        change_type: Type of change detected
        class_v1: Classification in version 1
        class_v2: Classification in version 2
        length_v1: Length in version 1 (0 if added)
        length_v2: Length in version 2 (0 if removed)
        size_change: Absolute size change in bp
        size_change_pct: Percent size change relative to v1
        notes: Additional notes about the change
    """

    name_v1: str | None
    name_v2: str | None
    change_type: ChangeType
    class_v1: str | None
    class_v2: str | None
    length_v1: int
    length_v2: int
    size_change: int = 0
    size_change_pct: float = 0.0
    notes: str = ""

    @property
    def primary_name(self) -> str:
        """Return the primary name (prefer v2, fall back to v1)."""
        return self.name_v2 or self.name_v1 or "unknown"


@dataclass
class VersionComparisonResult:
    """Results from comparing two assembly versions.

    Attributes:
        version1_label: Label for version 1 (e.g., filename or version number)
        version2_label: Label for version 2
        stats_v1: Assembly statistics for version 1
        stats_v2: Assembly statistics for version 2
        changes: List of detected scaffold changes
        unchanged_count: Number of scaffolds unchanged
        promoted_count: Number of scaffolds promoted
        demoted_count: Number of scaffolds demoted
        added_count: Number of new scaffolds
        removed_count: Number of removed scaffolds
        resized_count: Number of significantly resized scaffolds
        summary: Summary string of changes
    """

    version1_label: str
    version2_label: str
    stats_v1: AssemblyStats
    stats_v2: AssemblyStats
    changes: list[ScaffoldChange] = field(default_factory=list)
    unchanged_count: int = 0
    promoted_count: int = 0
    demoted_count: int = 0
    added_count: int = 0
    removed_count: int = 0
    resized_count: int = 0
    summary: str = ""

    @property
    def n50_change(self) -> int:
        """Change in N50 between versions."""
        return self.stats_v2.n50 - self.stats_v1.n50

    @property
    def n50_change_pct(self) -> float:
        """Percent change in N50."""
        if self.stats_v1.n50 == 0:
            return 0.0
        return (self.n50_change / self.stats_v1.n50) * 100

    @property
    def chromosome_count_change(self) -> int:
        """Change in chromosome count between versions."""
        return self.stats_v2.chromosome_count - self.stats_v1.chromosome_count

    @property
    def total_length_change(self) -> int:
        """Change in total assembly length."""
        return self.stats_v2.total_length - self.stats_v1.total_length

    @property
    def scaffold_count_change(self) -> int:
        """Change in total scaffold count."""
        return self.stats_v2.total_scaffolds - self.stats_v1.total_scaffolds


def _match_scaffolds(
    scaffolds_v1: list[ScaffoldInfo],
    scaffolds_v2: list[ScaffoldInfo],
    size_tolerance: float = 0.1,
) -> tuple[
    dict[str, str],  # name_v1 -> name_v2 for renamed scaffolds
    set[str],  # unmatched in v1 (removed)
    set[str],  # unmatched in v2 (added)
]:
    """
    Match scaffolds between two versions.

    First matches by exact name, then attempts to match unmatched scaffolds
    by similar size and chromosome ID.

    Args:
        scaffolds_v1: Scaffolds from version 1
        scaffolds_v2: Scaffolds from version 2
        size_tolerance: Maximum relative size difference for matching (0.1 = 10%)

    Returns:
        Tuple of (name_mapping, removed_names, added_names)
    """
    v1_by_name = {s.name: s for s in scaffolds_v1}
    v2_by_name = {s.name: s for s in scaffolds_v2}

    # Exact name matches
    name_mapping: dict[str, str] = {}
    matched_v1: set[str] = set()
    matched_v2: set[str] = set()

    for name in v1_by_name:
        if name in v2_by_name:
            name_mapping[name] = name
            matched_v1.add(name)
            matched_v2.add(name)

    # Unmatched scaffolds
    unmatched_v1 = [s for s in scaffolds_v1 if s.name not in matched_v1]
    unmatched_v2 = [s for s in scaffolds_v2 if s.name not in matched_v2]

    # Try to match by similar size and chromosome ID
    for s1 in unmatched_v1:
        best_match = None
        best_score = 0.0

        for s2 in unmatched_v2:
            if s2.name in matched_v2:
                continue

            # Calculate similarity score
            score = 0.0

            # Size similarity (within tolerance)
            if s1.length > 0:
                size_diff = abs(s1.length - s2.length) / s1.length
                if size_diff <= size_tolerance:
                    score += 1.0 - size_diff

            # Chromosome ID match is a strong signal
            if s1.chromosome_id and s1.chromosome_id == s2.chromosome_id:
                score += 2.0

            # Same classification type
            if s1.classification == s2.classification:
                score += 0.5

            if score > best_score and score >= 1.0:  # Minimum threshold
                best_score = score
                best_match = s2

        if best_match:
            name_mapping[s1.name] = best_match.name
            matched_v1.add(s1.name)
            matched_v2.add(best_match.name)

    removed = {s.name for s in scaffolds_v1 if s.name not in matched_v1}
    added = {s.name for s in scaffolds_v2 if s.name not in matched_v2}

    return name_mapping, removed, added


def _detect_change_type(
    s1: ScaffoldInfo | None,
    s2: ScaffoldInfo | None,
    is_renamed: bool = False,
    size_threshold: float = 0.2,
) -> tuple[ChangeType, str]:
    """
    Detect the type of change between two scaffold versions.

    Args:
        s1: Scaffold from version 1 (None if added)
        s2: Scaffold from version 2 (None if removed)
        is_renamed: Whether the scaffold was matched by heuristics (renamed)
        size_threshold: Minimum relative size change to flag as resized (0.2 = 20%)

    Returns:
        Tuple of (ChangeType, notes)
    """
    if s1 is None and s2 is not None:
        return ChangeType.ADDED, f"New scaffold in v2: {s2.classification}"

    if s1 is not None and s2 is None:
        return ChangeType.REMOVED, f"Removed from v2: was {s1.classification}"

    if s1 is None or s2 is None:
        return ChangeType.UNCHANGED, ""

    notes_parts = []

    # Check for significant size change
    size_change_pct = 0.0
    if s1.length > 0:
        size_change_pct = (s2.length - s1.length) / s1.length

    significant_resize = abs(size_change_pct) >= size_threshold

    # Check classification change
    rank1 = CLASSIFICATION_RANK.get(s1.classification, 0)
    rank2 = CLASSIFICATION_RANK.get(s2.classification, 0)

    if rank2 > rank1:
        # Promotion
        notes_parts.append(f"{s1.classification} → {s2.classification}")
        if significant_resize:
            if size_change_pct > 0:
                notes_parts.append(f"grew {size_change_pct*100:.1f}% (possible merge)")
            else:
                notes_parts.append(f"shrunk {abs(size_change_pct)*100:.1f}%")
        return ChangeType.PROMOTED, "; ".join(notes_parts)

    if rank2 < rank1:
        # Demotion
        notes_parts.append(f"{s1.classification} → {s2.classification}")
        if significant_resize:
            if size_change_pct < 0:
                notes_parts.append(
                    f"shrunk {abs(size_change_pct)*100:.1f}% (possible split)"
                )
            else:
                notes_parts.append(f"grew {size_change_pct*100:.1f}%")
        return ChangeType.DEMOTED, "; ".join(notes_parts)

    # Same rank - check for other changes
    if is_renamed:
        notes_parts.append(f"renamed: {s1.name} → {s2.name}")
        if significant_resize:
            notes_parts.append(f"size change: {size_change_pct*100:+.1f}%")
        return ChangeType.RENAMED, "; ".join(notes_parts)

    if significant_resize:
        if size_change_pct > 0:
            notes_parts.append(
                f"grew {size_change_pct*100:.1f}% ({s2.length - s1.length:,} bp)"
            )
        else:
            notes_parts.append(
                f"shrunk {abs(size_change_pct)*100:.1f}% ({s1.length - s2.length:,} bp)"
            )
        return ChangeType.RESIZED, "; ".join(notes_parts)

    if s1.classification != s2.classification:
        notes_parts.append(f"{s1.classification} → {s2.classification}")
        return ChangeType.RECLASSIFIED, "; ".join(notes_parts)

    return ChangeType.UNCHANGED, ""


def compare_assemblies(
    scaffolds_v1: list[ScaffoldInfo],
    stats_v1: AssemblyStats,
    scaffolds_v2: list[ScaffoldInfo],
    stats_v2: AssemblyStats,
    label_v1: str = "Version 1",
    label_v2: str = "Version 2",
    size_tolerance: float = 0.1,
    size_threshold: float = 0.2,
) -> VersionComparisonResult:
    """
    Compare two versions of an assembly.

    Args:
        scaffolds_v1: Classified scaffolds from version 1
        stats_v1: Assembly statistics for version 1
        scaffolds_v2: Classified scaffolds from version 2
        stats_v2: Assembly statistics for version 2
        label_v1: Label for version 1
        label_v2: Label for version 2
        size_tolerance: Max relative size difference for matching renamed scaffolds
        size_threshold: Min relative size change to flag as resized

    Returns:
        VersionComparisonResult with all detected changes
    """
    # Match scaffolds between versions
    name_mapping, removed_names, added_names = _match_scaffolds(
        scaffolds_v1, scaffolds_v2, size_tolerance
    )

    v1_by_name = {s.name: s for s in scaffolds_v1}
    v2_by_name = {s.name: s for s in scaffolds_v2}

    changes: list[ScaffoldChange] = []
    counts: dict[ChangeType, int] = dict.fromkeys(ChangeType, 0)

    # Process matched scaffolds
    for name_v1, name_v2 in name_mapping.items():
        s1 = v1_by_name[name_v1]
        s2 = v2_by_name[name_v2]
        is_renamed = name_v1 != name_v2

        change_type, notes = _detect_change_type(
            s1, s2, is_renamed=is_renamed, size_threshold=size_threshold
        )

        size_change = s2.length - s1.length
        size_change_pct = (size_change / s1.length * 100) if s1.length > 0 else 0.0

        changes.append(
            ScaffoldChange(
                name_v1=name_v1,
                name_v2=name_v2,
                change_type=change_type,
                class_v1=s1.classification,
                class_v2=s2.classification,
                length_v1=s1.length,
                length_v2=s2.length,
                size_change=size_change,
                size_change_pct=size_change_pct,
                notes=notes,
            )
        )
        counts[change_type] += 1

    # Process removed scaffolds
    for name in removed_names:
        s1 = v1_by_name[name]
        changes.append(
            ScaffoldChange(
                name_v1=name,
                name_v2=None,
                change_type=ChangeType.REMOVED,
                class_v1=s1.classification,
                class_v2=None,
                length_v1=s1.length,
                length_v2=0,
                size_change=-s1.length,
                size_change_pct=-100.0,
                notes=f"Removed from {label_v2}",
            )
        )
        counts[ChangeType.REMOVED] += 1

    # Process added scaffolds
    for name in added_names:
        s2 = v2_by_name[name]
        changes.append(
            ScaffoldChange(
                name_v1=None,
                name_v2=name,
                change_type=ChangeType.ADDED,
                class_v1=None,
                class_v2=s2.classification,
                length_v1=0,
                length_v2=s2.length,
                size_change=s2.length,
                size_change_pct=100.0,
                notes=f"Added in {label_v2}",
            )
        )
        counts[ChangeType.ADDED] += 1

    # Sort changes by significance (promotions/demotions first, then by size)
    priority_order = {
        ChangeType.PROMOTED: 0,
        ChangeType.DEMOTED: 1,
        ChangeType.RESIZED: 2,
        ChangeType.RENAMED: 3,
        ChangeType.ADDED: 4,
        ChangeType.REMOVED: 5,
        ChangeType.RECLASSIFIED: 6,
        ChangeType.UNCHANGED: 7,
    }
    changes.sort(key=lambda c: (priority_order[c.change_type], -abs(c.size_change)))

    # Generate summary
    summary_parts = []
    if counts[ChangeType.PROMOTED]:
        summary_parts.append(f"{counts[ChangeType.PROMOTED]} scaffolds promoted")
    if counts[ChangeType.DEMOTED]:
        summary_parts.append(f"{counts[ChangeType.DEMOTED]} scaffolds demoted")
    if counts[ChangeType.RESIZED]:
        summary_parts.append(f"{counts[ChangeType.RESIZED]} scaffolds resized")
    if counts[ChangeType.ADDED]:
        summary_parts.append(f"{counts[ChangeType.ADDED]} scaffolds added")
    if counts[ChangeType.REMOVED]:
        summary_parts.append(f"{counts[ChangeType.REMOVED]} scaffolds removed")

    n50_change = stats_v2.n50 - stats_v1.n50
    if n50_change != 0:
        n50_pct = (n50_change / stats_v1.n50 * 100) if stats_v1.n50 > 0 else 0
        summary_parts.append(f"N50: {n50_change:+,} bp ({n50_pct:+.1f}%)")

    chr_change = stats_v2.chromosome_count - stats_v1.chromosome_count
    if chr_change != 0:
        summary_parts.append(f"Chromosomes: {chr_change:+d}")

    summary = "; ".join(summary_parts) if summary_parts else "No significant changes"

    return VersionComparisonResult(
        version1_label=label_v1,
        version2_label=label_v2,
        stats_v1=stats_v1,
        stats_v2=stats_v2,
        changes=changes,
        unchanged_count=counts[ChangeType.UNCHANGED],
        promoted_count=counts[ChangeType.PROMOTED],
        demoted_count=counts[ChangeType.DEMOTED],
        added_count=counts[ChangeType.ADDED],
        removed_count=counts[ChangeType.REMOVED],
        resized_count=counts[ChangeType.RESIZED],
        summary=summary,
    )


def compare_fasta_files(
    fasta_v1: Path | str,
    fasta_v2: Path | str,
    label_v1: str | None = None,
    label_v2: str | None = None,
    min_chromosome_size: int = 10_000_000,
    size_tolerance: float = 0.1,
    size_threshold: float = 0.2,
) -> VersionComparisonResult:
    """
    Compare two FASTA files representing assembly versions.

    Args:
        fasta_v1: Path to version 1 FASTA file
        fasta_v2: Path to version 2 FASTA file
        label_v1: Label for version 1 (defaults to filename)
        label_v2: Label for version 2 (defaults to filename)
        min_chromosome_size: Minimum size for chromosome classification
        size_tolerance: Max relative size difference for matching
        size_threshold: Min relative size change to flag as resized

    Returns:
        VersionComparisonResult with all detected changes
    """
    fasta_v1 = Path(fasta_v1)
    fasta_v2 = Path(fasta_v2)

    # Use filenames as labels if not provided
    if label_v1 is None:
        label_v1 = fasta_v1.name
    if label_v2 is None:
        label_v2 = fasta_v2.name

    # Parse and classify both versions
    scaffolds_raw_v1 = parse_fasta(fasta_v1)
    scaffolds_raw_v2 = parse_fasta(fasta_v2)

    scaffolds_v1, stats_v1 = classify_scaffolds(
        scaffolds_raw_v1, min_chromosome_size=min_chromosome_size
    )
    scaffolds_v2, stats_v2 = classify_scaffolds(
        scaffolds_raw_v2, min_chromosome_size=min_chromosome_size
    )

    return compare_assemblies(
        scaffolds_v1,
        stats_v1,
        scaffolds_v2,
        stats_v2,
        label_v1=label_v1,
        label_v2=label_v2,
        size_tolerance=size_tolerance,
        size_threshold=size_threshold,
    )


def format_version_summary(result: VersionComparisonResult) -> str:
    """
    Format version comparison result as human-readable summary.

    Args:
        result: VersionComparisonResult from compare_assemblies()

    Returns:
        Formatted summary string
    """
    lines = []
    lines.append("=" * 70)
    lines.append("ASSEMBLY VERSION COMPARISON")
    lines.append("=" * 70)
    lines.append(f"Version 1: {result.version1_label}")
    lines.append(f"Version 2: {result.version2_label}")
    lines.append("")

    # Overall metrics comparison
    lines.append("METRIC CHANGES:")
    lines.append("-" * 40)
    lines.append(
        f"Total scaffolds:  {result.stats_v1.total_scaffolds:>10,} → "
        f"{result.stats_v2.total_scaffolds:>10,} "
        f"({result.scaffold_count_change:+d})"
    )
    lines.append(
        f"Total length:     {result.stats_v1.total_length:>10,} → "
        f"{result.stats_v2.total_length:>10,} "
        f"({result.total_length_change:+,} bp)"
    )
    lines.append(
        f"Chromosomes:      {result.stats_v1.chromosome_count:>10,} → "
        f"{result.stats_v2.chromosome_count:>10,} "
        f"({result.chromosome_count_change:+d})"
    )
    lines.append(
        f"N50:              {result.stats_v1.n50:>10,} → "
        f"{result.stats_v2.n50:>10,} "
        f"({result.n50_change:+,} bp, {result.n50_change_pct:+.1f}%)"
    )
    lines.append("")

    # Change summary
    lines.append("SCAFFOLD CHANGES:")
    lines.append("-" * 40)
    lines.append(f"Unchanged:   {result.unchanged_count:>5}")
    lines.append(f"Promoted:    {result.promoted_count:>5}")
    lines.append(f"Demoted:     {result.demoted_count:>5}")
    lines.append(f"Resized:     {result.resized_count:>5}")
    lines.append(f"Added:       {result.added_count:>5}")
    lines.append(f"Removed:     {result.removed_count:>5}")
    lines.append("")

    # Significant changes (exclude unchanged)
    significant = [c for c in result.changes if c.change_type != ChangeType.UNCHANGED]

    if significant:
        lines.append("SIGNIFICANT CHANGES:")
        lines.append("-" * 70)

        for change in significant[:50]:  # Limit to first 50
            change_icon = {
                ChangeType.PROMOTED: "↑",
                ChangeType.DEMOTED: "↓",
                ChangeType.RESIZED: "↔",
                ChangeType.RENAMED: "→",
                ChangeType.ADDED: "+",
                ChangeType.REMOVED: "-",
                ChangeType.RECLASSIFIED: "~",
            }.get(change.change_type, "?")

            name = change.primary_name
            lines.append(f"  {change_icon} {name}")
            if change.notes:
                lines.append(f"      {change.notes}")

        if len(significant) > 50:
            lines.append(f"  ... and {len(significant) - 50} more changes")

    lines.append("")
    lines.append(f"Summary: {result.summary}")
    lines.append("=" * 70)

    return "\n".join(lines)


def format_version_tsv(result: VersionComparisonResult) -> str:
    """
    Format version comparison result as TSV.

    Args:
        result: VersionComparisonResult from compare_assemblies()

    Returns:
        TSV-formatted string
    """
    lines = []

    # Header
    lines.append(
        "\t".join(
            [
                "change_type",
                "name_v1",
                "name_v2",
                "class_v1",
                "class_v2",
                "length_v1",
                "length_v2",
                "size_change",
                "size_change_pct",
                "notes",
            ]
        )
    )

    # Data rows
    for change in result.changes:
        lines.append(
            "\t".join(
                [
                    change.change_type.value,
                    change.name_v1 or "",
                    change.name_v2 or "",
                    change.class_v1 or "",
                    change.class_v2 or "",
                    str(change.length_v1),
                    str(change.length_v2),
                    str(change.size_change),
                    f"{change.size_change_pct:.2f}",
                    change.notes,
                ]
            )
        )

    return "\n".join(lines)


def format_version_json(result: VersionComparisonResult) -> dict:
    """
    Format version comparison result as JSON-serializable dict.

    Args:
        result: VersionComparisonResult from compare_assemblies()

    Returns:
        Dictionary suitable for JSON serialization
    """
    return {
        "version1": result.version1_label,
        "version2": result.version2_label,
        "stats_v1": result.stats_v1.to_dict(),
        "stats_v2": result.stats_v2.to_dict(),
        "metrics": {
            "n50_change": result.n50_change,
            "n50_change_pct": result.n50_change_pct,
            "chromosome_count_change": result.chromosome_count_change,
            "total_length_change": result.total_length_change,
            "scaffold_count_change": result.scaffold_count_change,
        },
        "counts": {
            "unchanged": result.unchanged_count,
            "promoted": result.promoted_count,
            "demoted": result.demoted_count,
            "resized": result.resized_count,
            "added": result.added_count,
            "removed": result.removed_count,
        },
        "changes": [
            {
                "name_v1": c.name_v1,
                "name_v2": c.name_v2,
                "change_type": c.change_type.value,
                "class_v1": c.class_v1,
                "class_v2": c.class_v2,
                "length_v1": c.length_v1,
                "length_v2": c.length_v2,
                "size_change": c.size_change,
                "size_change_pct": c.size_change_pct,
                "notes": c.notes,
            }
            for c in result.changes
        ],
        "summary": result.summary,
    }
