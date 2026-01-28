"""
Assembly Report Validation for ChromDetect.

This module validates that FASTA files match their corresponding NCBI
assembly reports, detecting discrepancies such as:
- Missing sequences (in report but not FASTA, or vice versa)
- Size mismatches between FASTA and report
- Name/accession matching issues
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from chromdetect.assembly_report import AssemblyReport, AssemblyReportEntry, parse_assembly_report
from chromdetect.core import parse_fasta


class ValidationSeverity(Enum):
    """Severity level for validation issues."""

    ERROR = "error"  # Critical issue that likely indicates a problem
    WARNING = "warning"  # Potential issue that should be reviewed
    INFO = "info"  # Informational note


class ValidationIssueType(Enum):
    """Types of validation issues that can be detected."""

    MISSING_IN_FASTA = "missing_in_fasta"  # Sequence in report but not in FASTA
    MISSING_IN_REPORT = "missing_in_report"  # Sequence in FASTA but not in report
    SIZE_MISMATCH = "size_mismatch"  # Sequence length differs between FASTA and report
    SIZE_NOT_IN_REPORT = "size_not_in_report"  # Report doesn't have size info


@dataclass
class ValidationIssue:
    """A single validation issue found during comparison.

    Attributes:
        issue_type: The type of issue detected
        severity: How severe the issue is
        sequence_name: Name of the affected sequence
        message: Human-readable description of the issue
        fasta_value: Value from FASTA (e.g., length) if applicable
        report_value: Value from report if applicable
    """

    issue_type: ValidationIssueType
    severity: ValidationSeverity
    sequence_name: str
    message: str
    fasta_value: int | str | None = None
    report_value: int | str | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "issue_type": self.issue_type.value,
            "severity": self.severity.value,
            "sequence_name": self.sequence_name,
            "message": self.message,
            "fasta_value": self.fasta_value,
            "report_value": self.report_value,
        }


@dataclass
class SequenceMatch:
    """A successfully matched sequence between FASTA and report.

    Attributes:
        fasta_name: Name in the FASTA file
        report_entry: Corresponding entry from assembly report
        match_type: How the match was made (name, genbank_accession, refseq_accession)
        fasta_length: Length from FASTA
    """

    fasta_name: str
    report_entry: AssemblyReportEntry
    match_type: str
    fasta_length: int

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "fasta_name": self.fasta_name,
            "report_name": self.report_entry.sequence_name,
            "match_type": self.match_type,
            "fasta_length": self.fasta_length,
            "report_length": self.report_entry.length,
            "sequence_role": self.report_entry.sequence_role,
            "assigned_molecule": self.report_entry.assigned_molecule,
        }


@dataclass
class ValidationResult:
    """Result of validating a FASTA against an assembly report.

    Attributes:
        fasta_path: Path to the FASTA file
        report_path: Path to the assembly report
        assembly_name: Name of the assembly from report
        organism: Organism name from report
        is_valid: True if no errors were found
        issues: List of validation issues found
        matches: List of successfully matched sequences
        fasta_only: Sequences in FASTA but not matched to report
        report_only: Sequences in report but not found in FASTA
        total_fasta_sequences: Total sequences in FASTA
        total_report_sequences: Total sequences in report
        matched_count: Number of sequences successfully matched
    """

    fasta_path: str
    report_path: str
    assembly_name: str | None
    organism: str | None
    is_valid: bool
    issues: list[ValidationIssue] = field(default_factory=list)
    matches: list[SequenceMatch] = field(default_factory=list)
    fasta_only: list[str] = field(default_factory=list)
    report_only: list[str] = field(default_factory=list)
    total_fasta_sequences: int = 0
    total_report_sequences: int = 0
    matched_count: int = 0

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "fasta_path": self.fasta_path,
            "report_path": self.report_path,
            "assembly_name": self.assembly_name,
            "organism": self.organism,
            "is_valid": self.is_valid,
            "summary": {
                "total_fasta_sequences": self.total_fasta_sequences,
                "total_report_sequences": self.total_report_sequences,
                "matched_count": self.matched_count,
                "fasta_only_count": len(self.fasta_only),
                "report_only_count": len(self.report_only),
                "error_count": len([i for i in self.issues if i.severity == ValidationSeverity.ERROR]),
                "warning_count": len([i for i in self.issues if i.severity == ValidationSeverity.WARNING]),
            },
            "issues": [i.to_dict() for i in self.issues],
            "matches": [m.to_dict() for m in self.matches],
            "fasta_only": self.fasta_only,
            "report_only": self.report_only,
        }

    @property
    def error_count(self) -> int:
        """Count of error-level issues."""
        return len([i for i in self.issues if i.severity == ValidationSeverity.ERROR])

    @property
    def warning_count(self) -> int:
        """Count of warning-level issues."""
        return len([i for i in self.issues if i.severity == ValidationSeverity.WARNING])


def _build_report_lookup(report: AssemblyReport) -> dict[str, AssemblyReportEntry]:
    """Build a lookup dictionary from all possible names/accessions to entries.

    Args:
        report: Parsed AssemblyReport

    Returns:
        Dictionary mapping sequence names and accessions to their entries
    """
    lookup: dict[str, AssemblyReportEntry] = {}

    for entry in report.entries:
        # Primary name
        lookup[entry.sequence_name] = entry

        # GenBank accession
        if entry.genbank_accession and entry.genbank_accession.lower() != "na":
            lookup[entry.genbank_accession] = entry

        # RefSeq accession
        if entry.refseq_accession and entry.refseq_accession.lower() != "na":
            lookup[entry.refseq_accession] = entry

    return lookup


def validate_fasta_against_report(
    fasta_path: Path | str,
    report_path: Path | str,
    size_tolerance: float = 0.0,
    strict: bool = False,
) -> ValidationResult:
    """
    Validate that a FASTA file matches an assembly report.

    Performs comprehensive validation including:
    - Checking all sequences in FASTA are present in report
    - Checking all sequences in report are present in FASTA
    - Verifying sequence lengths match (if reported)

    Args:
        fasta_path: Path to FASTA file
        report_path: Path to NCBI assembly report
        size_tolerance: Allowed fractional difference in size (0.0 = exact match)
        strict: If True, treat warnings as errors

    Returns:
        ValidationResult with detailed findings

    Raises:
        FileNotFoundError: If FASTA or report file doesn't exist
        ValueError: If files are invalid format
    """
    fasta_path = Path(fasta_path)
    report_path = Path(report_path)

    # Parse both files
    scaffolds = parse_fasta(fasta_path)
    report = parse_assembly_report(report_path)

    # Build lookups
    report_lookup = _build_report_lookup(report)

    # Track results
    issues: list[ValidationIssue] = []
    matches: list[SequenceMatch] = []
    matched_fasta_names: set[str] = set()
    matched_report_entries: set[str] = set()  # Track by sequence_name

    # Check each FASTA sequence against report
    for name, length, _seq in scaffolds:
        if name in report_lookup:
            entry = report_lookup[name]
            match_type = "name"

            # Determine match type
            if name == entry.sequence_name:
                match_type = "sequence_name"
            elif name == entry.genbank_accession:
                match_type = "genbank_accession"
            elif name == entry.refseq_accession:
                match_type = "refseq_accession"

            matches.append(SequenceMatch(
                fasta_name=name,
                report_entry=entry,
                match_type=match_type,
                fasta_length=length,
            ))
            matched_fasta_names.add(name)
            matched_report_entries.add(entry.sequence_name)

            # Check size match
            if entry.length is not None:
                if size_tolerance == 0.0:
                    if length != entry.length:
                        issues.append(ValidationIssue(
                            issue_type=ValidationIssueType.SIZE_MISMATCH,
                            severity=ValidationSeverity.ERROR,
                            sequence_name=name,
                            message=f"Size mismatch: FASTA has {length:,} bp, report has {entry.length:,} bp",
                            fasta_value=length,
                            report_value=entry.length,
                        ))
                else:
                    diff = abs(length - entry.length) / max(length, entry.length)
                    if diff > size_tolerance:
                        issues.append(ValidationIssue(
                            issue_type=ValidationIssueType.SIZE_MISMATCH,
                            severity=ValidationSeverity.ERROR,
                            sequence_name=name,
                            message=f"Size mismatch ({diff*100:.1f}% difference): FASTA has {length:,} bp, report has {entry.length:,} bp",
                            fasta_value=length,
                            report_value=entry.length,
                        ))
            else:
                issues.append(ValidationIssue(
                    issue_type=ValidationIssueType.SIZE_NOT_IN_REPORT,
                    severity=ValidationSeverity.INFO,
                    sequence_name=name,
                    message=f"Report does not include size for this sequence (FASTA: {length:,} bp)",
                    fasta_value=length,
                    report_value=None,
                ))
        else:
            # Not found in report
            issues.append(ValidationIssue(
                issue_type=ValidationIssueType.MISSING_IN_REPORT,
                severity=ValidationSeverity.WARNING if not strict else ValidationSeverity.ERROR,
                sequence_name=name,
                message=f"Sequence '{name}' in FASTA but not found in assembly report",
                fasta_value=length,
                report_value=None,
            ))

    # Check for sequences in report but not in FASTA
    for entry in report.entries:
        if entry.sequence_name not in matched_report_entries:
            # Check if we matched via accession
            matched_via_accession = False
            if entry.genbank_accession and entry.genbank_accession.lower() != "na":
                if entry.genbank_accession in matched_fasta_names:
                    matched_via_accession = True
            if entry.refseq_accession and entry.refseq_accession.lower() != "na":
                if entry.refseq_accession in matched_fasta_names:
                    matched_via_accession = True

            if not matched_via_accession:
                issues.append(ValidationIssue(
                    issue_type=ValidationIssueType.MISSING_IN_FASTA,
                    severity=ValidationSeverity.ERROR,
                    sequence_name=entry.sequence_name,
                    message=f"Sequence '{entry.sequence_name}' in report but not found in FASTA",
                    fasta_value=None,
                    report_value=entry.length,
                ))

    # Build result lists
    fasta_only = [name for name, _, _ in scaffolds if name not in matched_fasta_names]
    report_only = [
        entry.sequence_name for entry in report.entries
        if entry.sequence_name not in matched_report_entries
        and (entry.genbank_accession not in matched_fasta_names if entry.genbank_accession else True)
        and (entry.refseq_accession not in matched_fasta_names if entry.refseq_accession else True)
    ]

    # Determine overall validity
    error_count = len([i for i in issues if i.severity == ValidationSeverity.ERROR])
    is_valid = error_count == 0

    return ValidationResult(
        fasta_path=str(fasta_path),
        report_path=str(report_path),
        assembly_name=report.assembly_name,
        organism=report.organism,
        is_valid=is_valid,
        issues=issues,
        matches=matches,
        fasta_only=fasta_only,
        report_only=report_only,
        total_fasta_sequences=len(scaffolds),
        total_report_sequences=len(report.entries),
        matched_count=len(matches),
    )


def format_validation_summary(result: ValidationResult) -> str:
    """
    Format validation result as human-readable summary.

    Args:
        result: ValidationResult from validate_fasta_against_report

    Returns:
        Formatted string for display
    """
    lines = [
        "=" * 60,
        "ASSEMBLY REPORT VALIDATION",
        "=" * 60,
        "",
        f"FASTA file:      {result.fasta_path}",
        f"Assembly report: {result.report_path}",
    ]

    if result.assembly_name:
        lines.append(f"Assembly name:   {result.assembly_name}")
    if result.organism:
        lines.append(f"Organism:        {result.organism}")

    lines.extend([
        "",
        "-" * 60,
        "SUMMARY",
        "-" * 60,
        f"Total FASTA sequences:  {result.total_fasta_sequences:,}",
        f"Total report sequences: {result.total_report_sequences:,}",
        f"Successfully matched:   {result.matched_count:,}",
        f"In FASTA only:          {len(result.fasta_only):,}",
        f"In report only:         {len(result.report_only):,}",
        "",
    ])

    # Status
    if result.is_valid:
        lines.append("STATUS: VALID ✓")
    else:
        lines.append(f"STATUS: INVALID - {result.error_count} error(s) found")

    if result.warning_count > 0:
        lines.append(f"         {result.warning_count} warning(s)")

    # Issues
    if result.issues:
        lines.extend([
            "",
            "-" * 60,
            "ISSUES",
            "-" * 60,
        ])

        # Group by severity
        errors = [i for i in result.issues if i.severity == ValidationSeverity.ERROR]
        warnings = [i for i in result.issues if i.severity == ValidationSeverity.WARNING]
        infos = [i for i in result.issues if i.severity == ValidationSeverity.INFO]

        if errors:
            lines.append("")
            lines.append("ERRORS:")
            for issue in errors[:20]:  # Limit display
                lines.append(f"  [{issue.issue_type.value}] {issue.message}")
            if len(errors) > 20:
                lines.append(f"  ... and {len(errors) - 20} more errors")

        if warnings:
            lines.append("")
            lines.append("WARNINGS:")
            for issue in warnings[:20]:
                lines.append(f"  [{issue.issue_type.value}] {issue.message}")
            if len(warnings) > 20:
                lines.append(f"  ... and {len(warnings) - 20} more warnings")

        if infos and len(infos) <= 10:
            lines.append("")
            lines.append("INFO:")
            for issue in infos:
                lines.append(f"  [{issue.issue_type.value}] {issue.message}")
        elif infos:
            lines.append("")
            lines.append(f"INFO: {len(infos)} informational notes (use --format json for details)")

    return "\n".join(lines)


def format_validation_tsv(result: ValidationResult) -> str:
    """
    Format validation issues as TSV.

    Args:
        result: ValidationResult from validate_fasta_against_report

    Returns:
        TSV-formatted string
    """
    lines = ["severity\tissue_type\tsequence_name\tmessage\tfasta_value\treport_value"]

    for issue in result.issues:
        fasta_val = str(issue.fasta_value) if issue.fasta_value is not None else ""
        report_val = str(issue.report_value) if issue.report_value is not None else ""
        lines.append(
            f"{issue.severity.value}\t{issue.issue_type.value}\t"
            f"{issue.sequence_name}\t{issue.message}\t{fasta_val}\t{report_val}"
        )

    return "\n".join(lines)
