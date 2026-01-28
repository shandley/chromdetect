"""
Tests for assembly report validation functionality.
"""

import json

import pytest

from chromdetect.validation import (
    SequenceMatch,
    ValidationIssue,
    ValidationIssueType,
    ValidationResult,
    ValidationSeverity,
    format_validation_summary,
    format_validation_tsv,
    validate_fasta_against_report,
)


@pytest.fixture
def matching_fasta_content():
    """FASTA content that matches the sample report."""
    return """\
>chr1
ACGTACGTACGTACGT
>chr2
ACGTACGTACGTACGTACGT
>chrX
ACGTACGTACGTACGT
>chrM
ACGT
>chr1_unlocalized_001
ACGTACGT
>scaffold_001
ACGT
"""


@pytest.fixture
def matching_report_content():
    """Assembly report that matches the FASTA."""
    return """\
# Assembly name:  Test_Assembly_v1
# Organism name:  Test organism
# Taxid:          12345
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	16	chr1
chr2	assembled-molecule	2	Chromosome	CM000002.1	=	NC_000002.1	Primary Assembly	20	chr2
chrX	assembled-molecule	X	Chromosome	CM000023.1	=	NC_000023.1	Primary Assembly	16	chrX
chrM	assembled-molecule	MT	Mitochondrion	CM000024.1	=	NC_012920.1	non-nuclear	4	chrM
chr1_unlocalized_001	unlocalized-scaffold	1	Chromosome	CM000001.1	<>	na	Primary Assembly	8	na
scaffold_001	unplaced-scaffold	na	na	CM000025.1	=	na	Primary Assembly	4	na
"""


@pytest.fixture
def mismatched_size_fasta():
    """FASTA with size mismatches."""
    return """\
>chr1
ACGTACGTACGTACGTACGTACGT
>chr2
ACGT
"""


@pytest.fixture
def mismatched_size_report():
    """Report with different sizes."""
    return """\
# Assembly name:  Test_Assembly_v1
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	16	chr1
chr2	assembled-molecule	2	Chromosome	CM000002.1	=	NC_000002.1	Primary Assembly	100	chr2
"""


@pytest.fixture
def extra_fasta_sequences():
    """FASTA with sequences not in report."""
    return """\
>chr1
ACGTACGTACGTACGT
>extra_scaffold
ACGTACGT
>another_unknown
ACGT
"""


@pytest.fixture
def minimal_report():
    """Report with fewer sequences than FASTA."""
    return """\
# Assembly name:  Minimal_Assembly
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	16	chr1
"""


@pytest.fixture
def report_with_extra_sequences():
    """Report with sequences not in FASTA."""
    return """\
# Assembly name:  Test_Assembly
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	16	chr1
chr2	assembled-molecule	2	Chromosome	CM000002.1	=	NC_000002.1	Primary Assembly	20	chr2
chr3	assembled-molecule	3	Chromosome	CM000003.1	=	NC_000003.1	Primary Assembly	30	chr3
"""


@pytest.fixture
def accession_fasta():
    """FASTA using accession numbers instead of names."""
    return """\
>NC_000001.1
ACGTACGTACGTACGT
>CM000002.1
ACGTACGTACGTACGTACGT
"""


@pytest.fixture
def accession_report():
    """Report for accession matching."""
    return """\
# Assembly name:  Accession_Test
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	16	chr1
chr2	assembled-molecule	2	Chromosome	CM000002.1	=	NC_000002.1	Primary Assembly	20	chr2
"""


class TestValidationIssue:
    """Tests for ValidationIssue dataclass."""

    def test_issue_creation(self):
        """Test creating a validation issue."""
        issue = ValidationIssue(
            issue_type=ValidationIssueType.SIZE_MISMATCH,
            severity=ValidationSeverity.ERROR,
            sequence_name="chr1",
            message="Size mismatch",
            fasta_value=100,
            report_value=200,
        )
        assert issue.issue_type == ValidationIssueType.SIZE_MISMATCH
        assert issue.severity == ValidationSeverity.ERROR
        assert issue.fasta_value == 100

    def test_issue_to_dict(self):
        """Test serialization to dictionary."""
        issue = ValidationIssue(
            issue_type=ValidationIssueType.MISSING_IN_FASTA,
            severity=ValidationSeverity.ERROR,
            sequence_name="chr2",
            message="Missing sequence",
        )
        d = issue.to_dict()
        assert d["issue_type"] == "missing_in_fasta"
        assert d["severity"] == "error"
        assert d["sequence_name"] == "chr2"


class TestSequenceMatch:
    """Tests for SequenceMatch dataclass."""

    def test_match_to_dict(self):
        """Test serialization to dictionary."""
        from chromdetect.assembly_report import AssemblyReportEntry

        entry = AssemblyReportEntry(
            sequence_name="chr1",
            sequence_role="assembled-molecule",
            assigned_molecule="1",
            assigned_molecule_type="Chromosome",
            genbank_accession="CM000001.1",
            refseq_accession="NC_000001.1",
            length=100,
        )
        match = SequenceMatch(
            fasta_name="chr1",
            report_entry=entry,
            match_type="sequence_name",
            fasta_length=100,
        )
        d = match.to_dict()
        assert d["fasta_name"] == "chr1"
        assert d["report_name"] == "chr1"
        assert d["match_type"] == "sequence_name"
        assert d["fasta_length"] == 100


class TestValidationResult:
    """Tests for ValidationResult dataclass."""

    def test_result_to_dict(self):
        """Test serialization to dictionary."""
        result = ValidationResult(
            fasta_path="/path/to/file.fasta",
            report_path="/path/to/report.txt",
            assembly_name="Test",
            organism="Test organism",
            is_valid=True,
            total_fasta_sequences=10,
            total_report_sequences=10,
            matched_count=10,
        )
        d = result.to_dict()
        assert d["fasta_path"] == "/path/to/file.fasta"
        assert d["is_valid"] is True
        assert d["summary"]["total_fasta_sequences"] == 10

    def test_error_count_property(self):
        """Test error_count property."""
        result = ValidationResult(
            fasta_path="test.fasta",
            report_path="test.txt",
            assembly_name=None,
            organism=None,
            is_valid=False,
            issues=[
                ValidationIssue(
                    ValidationIssueType.SIZE_MISMATCH,
                    ValidationSeverity.ERROR,
                    "chr1",
                    "msg",
                ),
                ValidationIssue(
                    ValidationIssueType.MISSING_IN_REPORT,
                    ValidationSeverity.WARNING,
                    "chr2",
                    "msg",
                ),
                ValidationIssue(
                    ValidationIssueType.SIZE_MISMATCH,
                    ValidationSeverity.ERROR,
                    "chr3",
                    "msg",
                ),
            ],
        )
        assert result.error_count == 2
        assert result.warning_count == 1


class TestValidateFastaAgainstReport:
    """Tests for the main validation function."""

    def test_valid_matching_files(
        self, tmp_path, matching_fasta_content, matching_report_content
    ):
        """Test validation with perfectly matching files."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(matching_fasta_content)
        report_file.write_text(matching_report_content)

        result = validate_fasta_against_report(fasta_file, report_file)

        assert result.is_valid
        assert result.matched_count == 6
        assert result.error_count == 0
        assert len(result.fasta_only) == 0
        assert len(result.report_only) == 0

    def test_size_mismatch_detection(
        self, tmp_path, mismatched_size_fasta, mismatched_size_report
    ):
        """Test detection of size mismatches."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(mismatched_size_fasta)
        report_file.write_text(mismatched_size_report)

        result = validate_fasta_against_report(fasta_file, report_file)

        assert not result.is_valid
        size_issues = [
            i for i in result.issues
            if i.issue_type == ValidationIssueType.SIZE_MISMATCH
        ]
        assert len(size_issues) == 2  # Both chr1 and chr2 have wrong sizes

    def test_size_tolerance(
        self, tmp_path, mismatched_size_fasta, mismatched_size_report
    ):
        """Test that size tolerance allows small differences."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        # chr1: FASTA=24, report=16 (50% diff)
        # chr2: FASTA=4, report=100 (96% diff)
        fasta_file.write_text(mismatched_size_fasta)
        report_file.write_text(mismatched_size_report)

        # With 50% tolerance, chr1 should pass but chr2 should fail
        result = validate_fasta_against_report(
            fasta_file, report_file, size_tolerance=0.5
        )

        size_issues = [
            i for i in result.issues
            if i.issue_type == ValidationIssueType.SIZE_MISMATCH
        ]
        # chr2 has 96% difference, should still fail
        assert len(size_issues) == 1
        assert "chr2" in size_issues[0].sequence_name

    def test_missing_in_report_detection(
        self, tmp_path, extra_fasta_sequences, minimal_report
    ):
        """Test detection of sequences in FASTA but not in report."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(extra_fasta_sequences)
        report_file.write_text(minimal_report)

        result = validate_fasta_against_report(fasta_file, report_file)

        missing_issues = [
            i for i in result.issues
            if i.issue_type == ValidationIssueType.MISSING_IN_REPORT
        ]
        assert len(missing_issues) == 2  # extra_scaffold and another_unknown
        assert result.warning_count == 2  # These are warnings by default

    def test_missing_in_report_strict_mode(
        self, tmp_path, extra_fasta_sequences, minimal_report
    ):
        """Test that strict mode treats missing as errors."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(extra_fasta_sequences)
        report_file.write_text(minimal_report)

        result = validate_fasta_against_report(
            fasta_file, report_file, strict=True
        )

        assert not result.is_valid
        assert result.error_count == 2  # Now they're errors

    def test_missing_in_fasta_detection(
        self, tmp_path, report_with_extra_sequences
    ):
        """Test detection of sequences in report but not in FASTA."""
        fasta_content = """\
>chr1
ACGTACGTACGTACGT
"""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(fasta_content)
        report_file.write_text(report_with_extra_sequences)

        result = validate_fasta_against_report(fasta_file, report_file)

        assert not result.is_valid
        missing_issues = [
            i for i in result.issues
            if i.issue_type == ValidationIssueType.MISSING_IN_FASTA
        ]
        assert len(missing_issues) == 2  # chr2 and chr3

    def test_accession_matching(
        self, tmp_path, accession_fasta, accession_report
    ):
        """Test that accession numbers are matched correctly."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(accession_fasta)
        report_file.write_text(accession_report)

        result = validate_fasta_against_report(fasta_file, report_file)

        assert result.is_valid
        assert result.matched_count == 2

        # Check match types
        refseq_match = next(
            m for m in result.matches if m.fasta_name == "NC_000001.1"
        )
        assert refseq_match.match_type == "refseq_accession"

        genbank_match = next(
            m for m in result.matches if m.fasta_name == "CM000002.1"
        )
        assert genbank_match.match_type == "genbank_accession"

    def test_file_not_found_fasta(self, tmp_path):
        """Test error when FASTA file doesn't exist."""
        report_file = tmp_path / "report.txt"
        report_file.write_text("# test")

        with pytest.raises(FileNotFoundError):
            validate_fasta_against_report(
                tmp_path / "nonexistent.fasta", report_file
            )

    def test_file_not_found_report(self, tmp_path):
        """Test error when report file doesn't exist."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nACGT\n")

        with pytest.raises(FileNotFoundError):
            validate_fasta_against_report(
                fasta_file, tmp_path / "nonexistent.txt"
            )

    def test_report_without_sizes(self, tmp_path):
        """Test handling of report entries without sequence lengths."""
        fasta_content = ">chr1\nACGTACGT\n"
        report_content = """\
# Assembly name:  No_Sizes
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	na	chr1
"""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(fasta_content)
        report_file.write_text(report_content)

        result = validate_fasta_against_report(fasta_file, report_file)

        # Should be valid since we can't compare sizes
        assert result.is_valid
        # Should have an INFO issue about missing size
        info_issues = [
            i for i in result.issues
            if i.issue_type == ValidationIssueType.SIZE_NOT_IN_REPORT
        ]
        assert len(info_issues) == 1


class TestFormatValidationSummary:
    """Tests for summary formatting."""

    def test_valid_result_formatting(self, tmp_path, matching_fasta_content, matching_report_content):
        """Test formatting of a valid result."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(matching_fasta_content)
        report_file.write_text(matching_report_content)

        result = validate_fasta_against_report(fasta_file, report_file)
        summary = format_validation_summary(result)

        assert "VALID" in summary
        assert "Test_Assembly_v1" in summary
        assert "Test organism" in summary
        assert "Successfully matched:   6" in summary

    def test_invalid_result_formatting(self, tmp_path, mismatched_size_fasta, mismatched_size_report):
        """Test formatting of an invalid result."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(mismatched_size_fasta)
        report_file.write_text(mismatched_size_report)

        result = validate_fasta_against_report(fasta_file, report_file)
        summary = format_validation_summary(result)

        assert "INVALID" in summary
        assert "ERRORS:" in summary
        assert "size_mismatch" in summary


class TestFormatValidationTsv:
    """Tests for TSV formatting."""

    def test_tsv_output(self, tmp_path, mismatched_size_fasta, mismatched_size_report):
        """Test TSV format output."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(mismatched_size_fasta)
        report_file.write_text(mismatched_size_report)

        result = validate_fasta_against_report(fasta_file, report_file)
        tsv = format_validation_tsv(result)

        lines = tsv.strip().split("\n")
        assert lines[0] == "severity\tissue_type\tsequence_name\tmessage\tfasta_value\treport_value"
        assert len(lines) > 1  # Should have data rows


class TestValidationResultJson:
    """Tests for JSON serialization."""

    def test_json_serialization(self, tmp_path, matching_fasta_content, matching_report_content):
        """Test that result can be serialized to JSON."""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(matching_fasta_content)
        report_file.write_text(matching_report_content)

        result = validate_fasta_against_report(fasta_file, report_file)
        json_str = json.dumps(result.to_dict())

        # Should be valid JSON
        parsed = json.loads(json_str)
        assert parsed["is_valid"] is True
        assert parsed["summary"]["matched_count"] == 6


class TestEdgeCases:
    """Tests for edge cases and special scenarios."""

    def test_empty_fasta(self, tmp_path):
        """Test handling of empty FASTA file."""
        fasta_file = tmp_path / "empty.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text("")
        report_file.write_text("""\
# Assembly name:  Test
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	100	chr1
""")

        with pytest.raises(ValueError):
            validate_fasta_against_report(fasta_file, report_file)

    def test_gzipped_fasta(self, tmp_path, matching_report_content):
        """Test validation of gzipped FASTA files."""
        import gzip

        fasta_content = b">chr1\nACGTACGTACGTACGT\n"
        fasta_file = tmp_path / "test.fasta.gz"
        report_file = tmp_path / "report.txt"

        with gzip.open(fasta_file, "wb") as f:
            f.write(fasta_content)

        # Minimal report for chr1
        report_file.write_text("""\
# Assembly name:  Test
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	16	chr1
""")

        result = validate_fasta_against_report(fasta_file, report_file)
        assert result.is_valid
        assert result.matched_count == 1

    def test_special_characters_in_names(self, tmp_path):
        """Test handling of special characters in sequence names."""
        fasta_content = ">scaffold_1.1_random\nACGT\n"
        report_content = """\
# Assembly name:  Test
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
scaffold_1.1_random	unplaced-scaffold	na	na	CM000001.1	=	na	Primary Assembly	4	na
"""
        fasta_file = tmp_path / "test.fasta"
        report_file = tmp_path / "report.txt"
        fasta_file.write_text(fasta_content)
        report_file.write_text(report_content)

        result = validate_fasta_against_report(fasta_file, report_file)
        assert result.is_valid
        assert result.matched_count == 1
