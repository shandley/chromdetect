"""
Assembly Standardization for ChromDetect.

This module provides functionality to standardize scaffold/chromosome naming
conventions and convert between different naming systems used by various
databases and tools.

Supported naming conventions:
- UCSC: chr1, chr2, ..., chrX, chrY, chrM (human/mouse style)
- Ensembl: 1, 2, ..., X, Y, MT
- Simple: 1, 2, ..., X, Y, MT (same as Ensembl)
- GenBank: Accession-based (CM000001.1)
- RefSeq: Accession-based (NC_000001.11)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable

from chromdetect.core import parse_fasta, write_fasta


class NamingConvention(Enum):
    """Supported chromosome naming conventions."""

    UCSC = "ucsc"  # chr1, chr2, chrX, chrM
    ENSEMBL = "ensembl"  # 1, 2, X, MT
    SIMPLE = "simple"  # Same as Ensembl (1, 2, X, MT)
    GENBANK = "genbank"  # Accession-based (requires mapping)
    REFSEQ = "refseq"  # Accession-based (requires mapping)
    CUSTOM = "custom"  # User-defined mapping


# Standard chromosome ID normalization patterns
_CHR_PREFIX_PATTERN = re.compile(r"^(?:chr(?:omosome)?[_-]?)?(.+)$", re.IGNORECASE)

# Mitochondria aliases
_MT_ALIASES = {"m", "mt", "mito", "mitochondria", "mitochondrion", "chrm", "chrmt"}


def normalize_chromosome_id(name: str) -> str:
    """
    Normalize a chromosome name to its core identifier.

    Removes common prefixes like 'chr', 'chromosome', etc.

    Args:
        name: Chromosome name (e.g., 'chr1', 'chromosome_X', 'Chr_MT')

    Returns:
        Normalized ID (e.g., '1', 'X', 'MT')

    Examples:
        >>> normalize_chromosome_id('chr1')
        '1'
        >>> normalize_chromosome_id('chromosome_X')
        'X'
        >>> normalize_chromosome_id('chrM')
        'MT'
    """
    # Remove prefix
    match = _CHR_PREFIX_PATTERN.match(name)
    if match:
        core = match.group(1)
    else:
        core = name

    # Normalize case
    core_lower = core.lower()

    # Normalize mitochondria
    if core_lower in _MT_ALIASES:
        return "MT"

    # Keep numeric chromosomes as-is, uppercase letters
    if core.isdigit():
        return core
    elif len(core) <= 2 and core.isalpha():
        return core.upper()
    else:
        return core


def to_ucsc(chromosome_id: str) -> str:
    """
    Convert chromosome ID to UCSC naming convention.

    Args:
        chromosome_id: Normalized chromosome ID (e.g., '1', 'X', 'MT')

    Returns:
        UCSC-style name (e.g., 'chr1', 'chrX', 'chrM')
    """
    normalized = normalize_chromosome_id(chromosome_id)

    # MT -> M for UCSC
    if normalized == "MT":
        return "chrM"

    return f"chr{normalized}"


def to_ensembl(chromosome_id: str) -> str:
    """
    Convert chromosome ID to Ensembl naming convention.

    Args:
        chromosome_id: Chromosome ID (may have prefix)

    Returns:
        Ensembl-style name (e.g., '1', 'X', 'MT')
    """
    return normalize_chromosome_id(chromosome_id)


def to_simple(chromosome_id: str) -> str:
    """
    Convert chromosome ID to simple naming convention.

    Same as Ensembl style.

    Args:
        chromosome_id: Chromosome ID (may have prefix)

    Returns:
        Simple name (e.g., '1', 'X', 'MT')
    """
    return normalize_chromosome_id(chromosome_id)


# Human chromosome accession mappings (GRCh38)
HUMAN_REFSEQ_MAP: dict[str, str] = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}

HUMAN_GENBANK_MAP: dict[str, str] = {
    "1": "CM000663.2", "2": "CM000664.2", "3": "CM000665.2",
    "4": "CM000666.2", "5": "CM000667.2", "6": "CM000668.2",
    "7": "CM000669.2", "8": "CM000670.2", "9": "CM000671.2",
    "10": "CM000672.2", "11": "CM000673.2", "12": "CM000674.2",
    "13": "CM000675.2", "14": "CM000676.2", "15": "CM000677.2",
    "16": "CM000678.2", "17": "CM000679.2", "18": "CM000680.2",
    "19": "CM000681.2", "20": "CM000682.2", "21": "CM000683.2",
    "22": "CM000684.2", "X": "CM000685.2", "Y": "CM000686.2",
    "MT": "J01415.2",
}

# Mouse chromosome accession mappings (GRCm39)
MOUSE_REFSEQ_MAP: dict[str, str] = {
    "1": "NC_000067.7", "2": "NC_000068.8", "3": "NC_000069.7",
    "4": "NC_000070.7", "5": "NC_000071.7", "6": "NC_000072.7",
    "7": "NC_000073.7", "8": "NC_000074.7", "9": "NC_000075.7",
    "10": "NC_000076.7", "11": "NC_000077.7", "12": "NC_000078.7",
    "13": "NC_000079.7", "14": "NC_000080.7", "15": "NC_000081.7",
    "16": "NC_000082.7", "17": "NC_000083.7", "18": "NC_000084.7",
    "19": "NC_000085.7", "X": "NC_000086.8", "Y": "NC_000087.8",
    "MT": "NC_005089.1",
}

# Reverse mappings (accession -> chromosome ID)
HUMAN_REFSEQ_REVERSE = {v: k for k, v in HUMAN_REFSEQ_MAP.items()}
HUMAN_GENBANK_REVERSE = {v: k for k, v in HUMAN_GENBANK_MAP.items()}
MOUSE_REFSEQ_REVERSE = {v: k for k, v in MOUSE_REFSEQ_MAP.items()}


@dataclass
class ConversionMapping:
    """Mapping configuration for chromosome name conversion.

    Attributes:
        source_convention: Source naming convention
        target_convention: Target naming convention
        species: Species for accession-based conversions
        custom_map: Custom mapping dictionary (name -> new_name)
    """

    source_convention: NamingConvention | None = None
    target_convention: NamingConvention = NamingConvention.UCSC
    species: str = "human"
    custom_map: dict[str, str] = field(default_factory=dict)


def get_converter(
    target: NamingConvention | str,
    species: str = "human",
    custom_map: dict[str, str] | None = None,
) -> Callable[[str], str]:
    """
    Get a converter function for the specified target convention.

    Args:
        target: Target naming convention
        species: Species for accession-based conversions
        custom_map: Optional custom mapping dictionary

    Returns:
        Converter function that takes a chromosome name and returns converted name
    """
    if isinstance(target, str):
        try:
            target = NamingConvention(target.lower())
        except ValueError:
            raise ValueError(f"Unknown naming convention: {target}") from None

    if custom_map:
        def custom_converter(name: str) -> str:
            if name in custom_map:
                return custom_map[name]
            normalized = normalize_chromosome_id(name)
            if normalized in custom_map:
                return custom_map[normalized]
            return name
        return custom_converter

    if target == NamingConvention.UCSC:
        return to_ucsc
    elif target in (NamingConvention.ENSEMBL, NamingConvention.SIMPLE):
        return to_ensembl
    elif target == NamingConvention.REFSEQ:
        species_lower = species.lower()
        if species_lower in ("human", "homo sapiens", "hsapiens"):
            refseq_map = HUMAN_REFSEQ_MAP
        elif species_lower in ("mouse", "mus musculus", "mmusculus"):
            refseq_map = MOUSE_REFSEQ_MAP
        else:
            raise ValueError(f"RefSeq mapping not available for species: {species}")

        def refseq_converter(name: str) -> str:
            normalized = normalize_chromosome_id(name)
            return refseq_map.get(normalized, name)
        return refseq_converter

    elif target == NamingConvention.GENBANK:
        species_lower = species.lower()
        if species_lower in ("human", "homo sapiens", "hsapiens"):
            genbank_map = HUMAN_GENBANK_MAP
        else:
            raise ValueError(f"GenBank mapping not available for species: {species}")

        def genbank_converter(name: str) -> str:
            normalized = normalize_chromosome_id(name)
            return genbank_map.get(normalized, name)
        return genbank_converter

    else:
        raise ValueError(f"Unsupported target convention: {target}")


def detect_convention(names: list[str]) -> NamingConvention | None:
    """
    Attempt to detect the naming convention used in a list of chromosome names.

    Args:
        names: List of chromosome/scaffold names

    Returns:
        Detected NamingConvention or None if unable to determine
    """
    if not names:
        return None

    chr_prefix_count = 0
    ensembl_count = 0
    refseq_count = 0
    genbank_count = 0

    for name in names:
        name_lower = name.lower()

        if name_lower.startswith("chr"):
            chr_prefix_count += 1
        elif name.startswith("NC_"):
            refseq_count += 1
        elif name.startswith("CM"):
            genbank_count += 1
        elif re.match(r"^\d+$", name) or name.upper() in ("X", "Y", "MT", "W", "Z"):
            ensembl_count += 1

    total = len(names)
    threshold = 0.5  # At least 50% should match

    if refseq_count / total >= threshold:
        return NamingConvention.REFSEQ
    elif genbank_count / total >= threshold:
        return NamingConvention.GENBANK
    elif chr_prefix_count / total >= threshold:
        return NamingConvention.UCSC
    elif ensembl_count / total >= threshold:
        return NamingConvention.ENSEMBL

    return None


@dataclass
class RenameResult:
    """Result of renaming scaffolds in a FASTA file.

    Attributes:
        original_names: List of original scaffold names
        new_names: List of new scaffold names
        mapping: Dictionary of original -> new name
        unchanged: Names that weren't changed
        renamed_count: Number of names that were changed
    """

    original_names: list[str]
    new_names: list[str]
    mapping: dict[str, str]
    unchanged: list[str]
    renamed_count: int

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "original_names": self.original_names,
            "new_names": self.new_names,
            "mapping": self.mapping,
            "unchanged": self.unchanged,
            "renamed_count": self.renamed_count,
            "total_scaffolds": len(self.original_names),
        }


def rename_scaffolds(
    scaffolds: list[tuple[str, int, str]],
    target_convention: NamingConvention | str,
    species: str = "human",
    custom_map: dict[str, str] | None = None,
    chromosomes_only: bool = False,
) -> tuple[list[tuple[str, int, str]], RenameResult]:
    """
    Rename scaffolds according to a target naming convention.

    Args:
        scaffolds: List of (name, length, sequence) tuples from parse_fasta()
        target_convention: Target naming convention
        species: Species for accession-based conversions
        custom_map: Optional custom mapping dictionary
        chromosomes_only: If True, only rename chromosome-like scaffolds

    Returns:
        Tuple of (renamed_scaffolds, RenameResult)
    """
    converter = get_converter(target_convention, species, custom_map)

    original_names = []
    new_names = []
    mapping = {}
    unchanged = []
    renamed_scaffolds = []

    for name, length, seq in scaffolds:
        original_names.append(name)

        # Check if this looks like a chromosome
        normalized = normalize_chromosome_id(name)
        is_chromosome_like = (
            normalized.isdigit() or
            normalized.upper() in ("X", "Y", "Z", "W", "MT", "M") or
            len(normalized) <= 3
        )

        if chromosomes_only and not is_chromosome_like:
            new_name = name
            unchanged.append(name)
        else:
            new_name = converter(name)
            if new_name == name:
                unchanged.append(name)

        new_names.append(new_name)
        mapping[name] = new_name
        renamed_scaffolds.append((new_name, length, seq))

    renamed_count = sum(1 for old, new in zip(original_names, new_names) if old != new)

    result = RenameResult(
        original_names=original_names,
        new_names=new_names,
        mapping=mapping,
        unchanged=unchanged,
        renamed_count=renamed_count,
    )

    return renamed_scaffolds, result


def standardize_fasta(
    input_path: Path | str,
    output_path: Path | str,
    target_convention: NamingConvention | str,
    species: str = "human",
    custom_map: dict[str, str] | None = None,
    chromosomes_only: bool = False,
) -> RenameResult:
    """
    Read a FASTA file, rename scaffolds, and write to a new file.

    Args:
        input_path: Path to input FASTA file
        output_path: Path to output FASTA file
        target_convention: Target naming convention
        species: Species for accession-based conversions
        custom_map: Optional custom mapping dictionary
        chromosomes_only: If True, only rename chromosome-like scaffolds

    Returns:
        RenameResult with details of the conversion
    """
    scaffolds = parse_fasta(input_path, keep_full_sequence=True)
    renamed_scaffolds, result = rename_scaffolds(
        scaffolds,
        target_convention,
        species,
        custom_map,
        chromosomes_only,
    )

    sequences = [(name, seq) for name, _length, seq in renamed_scaffolds]
    write_fasta(sequences, output_path)

    return result


# NCBI submission compliance checking

@dataclass
class ComplianceIssue:
    """An issue found during NCBI compliance checking."""

    severity: str  # "error", "warning"
    issue_type: str
    scaffold_name: str
    message: str

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "severity": self.severity,
            "issue_type": self.issue_type,
            "scaffold_name": self.scaffold_name,
            "message": self.message,
        }


@dataclass
class ComplianceResult:
    """Result of NCBI compliance checking.

    Attributes:
        is_compliant: True if no errors found
        issues: List of compliance issues
        total_scaffolds: Total number of scaffolds checked
    """

    is_compliant: bool
    issues: list[ComplianceIssue] = field(default_factory=list)
    total_scaffolds: int = 0

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "is_compliant": self.is_compliant,
            "total_scaffolds": self.total_scaffolds,
            "error_count": len([i for i in self.issues if i.severity == "error"]),
            "warning_count": len([i for i in self.issues if i.severity == "warning"]),
            "issues": [i.to_dict() for i in self.issues],
        }


# NCBI forbidden characters in sequence IDs
NCBI_FORBIDDEN_CHARS = set(' \t\n\r>[]|"\'\\')
NCBI_MAX_ID_LENGTH = 50  # NCBI recommends IDs <= 50 characters


def check_ncbi_compliance(
    scaffolds: list[tuple[str, int, str]] | list[str],
) -> ComplianceResult:
    """
    Check if scaffold names comply with NCBI submission requirements.

    NCBI requirements:
    - No spaces or tabs in sequence IDs
    - No special characters: > [ ] | " ' \\
    - Unique names (no duplicates)
    - Recommended: IDs <= 50 characters
    - Recommended: Alphanumeric with underscores

    Args:
        scaffolds: List of (name, length, seq) tuples or just names

    Returns:
        ComplianceResult with findings
    """
    issues: list[ComplianceIssue] = []

    # Handle both full scaffold tuples and just names
    if scaffolds and isinstance(scaffolds[0], tuple):
        names = [s[0] for s in scaffolds]  # type: ignore
    else:
        names = scaffolds  # type: ignore

    seen_names: dict[str, int] = {}

    for i, name in enumerate(names):
        # Check for forbidden characters
        forbidden_found = [c for c in name if c in NCBI_FORBIDDEN_CHARS]
        if forbidden_found:
            issues.append(ComplianceIssue(
                severity="error",
                issue_type="forbidden_characters",
                scaffold_name=name,
                message=f"Contains forbidden characters: {forbidden_found!r}",
            ))

        # Check for empty name
        if not name or name.isspace():
            issues.append(ComplianceIssue(
                severity="error",
                issue_type="empty_name",
                scaffold_name=name or "(empty)",
                message="Scaffold name is empty or whitespace",
            ))
            continue

        # Check length
        if len(name) > NCBI_MAX_ID_LENGTH:
            issues.append(ComplianceIssue(
                severity="warning",
                issue_type="long_name",
                scaffold_name=name,
                message=f"Name exceeds recommended {NCBI_MAX_ID_LENGTH} characters ({len(name)} chars)",
            ))

        # Check for duplicates
        if name in seen_names:
            issues.append(ComplianceIssue(
                severity="error",
                issue_type="duplicate_name",
                scaffold_name=name,
                message=f"Duplicate name (first seen at position {seen_names[name] + 1})",
            ))
        else:
            seen_names[name] = i

        # Check for non-recommended characters (warning only)
        if not re.match(r"^[A-Za-z0-9_.-]+$", name):
            issues.append(ComplianceIssue(
                severity="warning",
                issue_type="non_standard_characters",
                scaffold_name=name,
                message="Contains non-standard characters (recommend: alphanumeric, underscore, dash, dot)",
            ))

    is_compliant = not any(i.severity == "error" for i in issues)

    return ComplianceResult(
        is_compliant=is_compliant,
        issues=issues,
        total_scaffolds=len(names),
    )


def format_compliance_summary(result: ComplianceResult) -> str:
    """
    Format compliance result as human-readable summary.

    Args:
        result: ComplianceResult from check_ncbi_compliance

    Returns:
        Formatted string for display
    """
    lines = [
        "=" * 60,
        "NCBI SUBMISSION COMPLIANCE CHECK",
        "=" * 60,
        "",
        f"Total scaffolds checked: {result.total_scaffolds}",
        "",
    ]

    error_count = len([i for i in result.issues if i.severity == "error"])
    warning_count = len([i for i in result.issues if i.severity == "warning"])

    if result.is_compliant:
        lines.append("STATUS: COMPLIANT")
    else:
        lines.append(f"STATUS: NOT COMPLIANT - {error_count} error(s)")

    if warning_count > 0:
        lines.append(f"         {warning_count} warning(s)")

    if result.issues:
        lines.extend([
            "",
            "-" * 60,
            "ISSUES",
            "-" * 60,
        ])

        errors = [i for i in result.issues if i.severity == "error"]
        warnings = [i for i in result.issues if i.severity == "warning"]

        if errors:
            lines.append("")
            lines.append("ERRORS:")
            for issue in errors[:20]:
                lines.append(f"  [{issue.issue_type}] {issue.scaffold_name}: {issue.message}")
            if len(errors) > 20:
                lines.append(f"  ... and {len(errors) - 20} more errors")

        if warnings:
            lines.append("")
            lines.append("WARNINGS:")
            for issue in warnings[:20]:
                lines.append(f"  [{issue.issue_type}] {issue.scaffold_name}: {issue.message}")
            if len(warnings) > 20:
                lines.append(f"  ... and {len(warnings) - 20} more warnings")

    return "\n".join(lines)


def format_rename_summary(result: RenameResult) -> str:
    """
    Format rename result as human-readable summary.

    Args:
        result: RenameResult from rename_scaffolds

    Returns:
        Formatted string for display
    """
    lines = [
        "=" * 60,
        "SCAFFOLD RENAMING SUMMARY",
        "=" * 60,
        "",
        f"Total scaffolds:  {len(result.original_names)}",
        f"Renamed:          {result.renamed_count}",
        f"Unchanged:        {len(result.unchanged)}",
        "",
    ]

    if result.renamed_count > 0:
        lines.extend([
            "-" * 60,
            "RENAMED SCAFFOLDS",
            "-" * 60,
        ])

        renamed = [(old, new) for old, new in result.mapping.items() if old != new]
        for old, new in renamed[:30]:
            lines.append(f"  {old} -> {new}")
        if len(renamed) > 30:
            lines.append(f"  ... and {len(renamed) - 30} more")

    return "\n".join(lines)
