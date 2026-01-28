"""
Tests for assembly standardization functionality.
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from chromdetect.standardize import (
    ComplianceResult,
    NamingConvention,
    RenameResult,
    check_ncbi_compliance,
    detect_convention,
    format_compliance_summary,
    format_rename_summary,
    get_converter,
    normalize_chromosome_id,
    rename_scaffolds,
    standardize_fasta,
    to_ensembl,
    to_ucsc,
)


class TestNormalizeChromosomeId:
    """Tests for chromosome ID normalization."""

    def test_chr_prefix(self):
        """Test removing chr prefix."""
        assert normalize_chromosome_id("chr1") == "1"
        assert normalize_chromosome_id("chr22") == "22"
        assert normalize_chromosome_id("chrX") == "X"
        assert normalize_chromosome_id("chrY") == "Y"

    def test_chromosome_prefix(self):
        """Test removing chromosome prefix."""
        assert normalize_chromosome_id("chromosome_1") == "1"
        assert normalize_chromosome_id("chromosome1") == "1"
        assert normalize_chromosome_id("Chromosome_X") == "X"

    def test_mitochondria_normalization(self):
        """Test mitochondria aliases normalize to MT."""
        assert normalize_chromosome_id("chrM") == "MT"
        assert normalize_chromosome_id("chrMT") == "MT"
        assert normalize_chromosome_id("M") == "MT"
        assert normalize_chromosome_id("MT") == "MT"
        assert normalize_chromosome_id("mito") == "MT"

    def test_already_normalized(self):
        """Test already normalized IDs."""
        assert normalize_chromosome_id("1") == "1"
        assert normalize_chromosome_id("22") == "22"
        assert normalize_chromosome_id("X") == "X"
        assert normalize_chromosome_id("MT") == "MT"

    def test_case_normalization(self):
        """Test case normalization for letters."""
        assert normalize_chromosome_id("x") == "X"
        assert normalize_chromosome_id("y") == "Y"
        assert normalize_chromosome_id("Chr_x") == "X"


class TestToUcsc:
    """Tests for UCSC conversion."""

    def test_basic_conversion(self):
        """Test basic UCSC conversion."""
        assert to_ucsc("1") == "chr1"
        assert to_ucsc("22") == "chr22"
        assert to_ucsc("X") == "chrX"
        assert to_ucsc("Y") == "chrY"

    def test_mt_conversion(self):
        """Test MT -> chrM conversion."""
        assert to_ucsc("MT") == "chrM"
        assert to_ucsc("chrMT") == "chrM"

    def test_already_ucsc(self):
        """Test already UCSC names."""
        assert to_ucsc("chr1") == "chr1"
        assert to_ucsc("chrX") == "chrX"


class TestToEnsembl:
    """Tests for Ensembl conversion."""

    def test_basic_conversion(self):
        """Test basic Ensembl conversion."""
        assert to_ensembl("chr1") == "1"
        assert to_ensembl("chr22") == "22"
        assert to_ensembl("chrX") == "X"

    def test_mt_conversion(self):
        """Test chrM -> MT conversion."""
        assert to_ensembl("chrM") == "MT"

    def test_already_ensembl(self):
        """Test already Ensembl names."""
        assert to_ensembl("1") == "1"
        assert to_ensembl("X") == "X"
        assert to_ensembl("MT") == "MT"


class TestGetConverter:
    """Tests for get_converter function."""

    def test_ucsc_converter(self):
        """Test getting UCSC converter."""
        converter = get_converter(NamingConvention.UCSC)
        assert converter("1") == "chr1"
        assert converter("X") == "chrX"

    def test_ensembl_converter(self):
        """Test getting Ensembl converter."""
        converter = get_converter(NamingConvention.ENSEMBL)
        assert converter("chr1") == "1"
        assert converter("chrX") == "X"

    def test_string_convention(self):
        """Test using string for convention."""
        converter = get_converter("ucsc")
        assert converter("1") == "chr1"

    def test_refseq_converter_human(self):
        """Test RefSeq converter for human."""
        converter = get_converter(NamingConvention.REFSEQ, species="human")
        assert converter("1") == "NC_000001.11"
        assert converter("X") == "NC_000023.11"
        assert converter("chr1") == "NC_000001.11"

    def test_refseq_converter_mouse(self):
        """Test RefSeq converter for mouse."""
        converter = get_converter(NamingConvention.REFSEQ, species="mouse")
        assert converter("1") == "NC_000067.7"
        assert converter("X") == "NC_000086.8"

    def test_genbank_converter_human(self):
        """Test GenBank converter for human."""
        converter = get_converter(NamingConvention.GENBANK, species="human")
        assert converter("1") == "CM000663.2"
        assert converter("X") == "CM000685.2"

    def test_custom_map(self):
        """Test custom mapping."""
        custom = {"chr1": "custom_chr1", "chr2": "custom_chr2"}
        converter = get_converter(NamingConvention.UCSC, custom_map=custom)
        assert converter("chr1") == "custom_chr1"
        assert converter("chr2") == "custom_chr2"
        assert converter("chr3") == "chr3"  # Not in custom map

    def test_invalid_species_refseq(self):
        """Test error for unsupported species."""
        with pytest.raises(ValueError, match="not available"):
            get_converter(NamingConvention.REFSEQ, species="unknown_species")

    def test_invalid_convention(self):
        """Test error for invalid convention."""
        with pytest.raises(ValueError, match="Unknown naming convention"):
            get_converter("invalid_convention")


class TestDetectConvention:
    """Tests for convention detection."""

    def test_detect_ucsc(self):
        """Test detecting UCSC convention."""
        names = ["chr1", "chr2", "chr3", "chrX", "chrY"]
        assert detect_convention(names) == NamingConvention.UCSC

    def test_detect_ensembl(self):
        """Test detecting Ensembl convention."""
        names = ["1", "2", "3", "X", "Y", "MT"]
        assert detect_convention(names) == NamingConvention.ENSEMBL

    def test_detect_refseq(self):
        """Test detecting RefSeq convention."""
        names = ["NC_000001.11", "NC_000002.12", "NC_000003.12"]
        assert detect_convention(names) == NamingConvention.REFSEQ

    def test_detect_genbank(self):
        """Test detecting GenBank convention."""
        names = ["CM000663.2", "CM000664.2", "CM000665.2"]
        assert detect_convention(names) == NamingConvention.GENBANK

    def test_empty_list(self):
        """Test empty list returns None."""
        assert detect_convention([]) is None

    def test_mixed_conventions(self):
        """Test mixed conventions may return None or most common."""
        names = ["chr1", "2", "NC_000003.12", "scaffold_1"]
        # Should return most common or None
        result = detect_convention(names)
        # No clear majority, so might be None
        assert result is None or isinstance(result, NamingConvention)


class TestRenameScaffolds:
    """Tests for scaffold renaming."""

    def test_rename_to_ucsc(self):
        """Test renaming to UCSC convention."""
        scaffolds = [
            ("1", 1000, "ACGT"),
            ("2", 2000, "GCTA"),
            ("X", 3000, "TTTT"),
        ]
        renamed, result = rename_scaffolds(scaffolds, NamingConvention.UCSC)

        assert renamed[0][0] == "chr1"
        assert renamed[1][0] == "chr2"
        assert renamed[2][0] == "chrX"
        assert result.renamed_count == 3

    def test_rename_to_ensembl(self):
        """Test renaming to Ensembl convention."""
        scaffolds = [
            ("chr1", 1000, "ACGT"),
            ("chr2", 2000, "GCTA"),
            ("chrM", 3000, "TTTT"),
        ]
        renamed, result = rename_scaffolds(scaffolds, NamingConvention.ENSEMBL)

        assert renamed[0][0] == "1"
        assert renamed[1][0] == "2"
        assert renamed[2][0] == "MT"
        assert result.renamed_count == 3

    def test_chromosomes_only(self):
        """Test renaming only chromosome-like scaffolds."""
        scaffolds = [
            ("chr1", 1000, "ACGT"),
            ("scaffold_123", 2000, "GCTA"),
            ("contig_abc", 3000, "TTTT"),
        ]
        renamed, result = rename_scaffolds(
            scaffolds,
            NamingConvention.ENSEMBL,
            chromosomes_only=True,
        )

        assert renamed[0][0] == "1"
        assert renamed[1][0] == "scaffold_123"  # Unchanged
        assert renamed[2][0] == "contig_abc"  # Unchanged
        assert result.renamed_count == 1

    def test_rename_result_structure(self):
        """Test RenameResult structure."""
        scaffolds = [("chr1", 1000, "ACGT")]
        _, result = rename_scaffolds(scaffolds, NamingConvention.ENSEMBL)

        assert "chr1" in result.original_names
        assert "1" in result.new_names
        assert result.mapping["chr1"] == "1"


class TestStandardizeFasta:
    """Tests for FASTA file standardization."""

    def test_standardize_file(self, tmp_path):
        """Test standardizing a FASTA file."""
        input_file = tmp_path / "input.fasta"
        output_file = tmp_path / "output.fasta"

        input_file.write_text(">1\nACGT\n>2\nGCTA\n>X\nTTTT\n")

        result = standardize_fasta(input_file, output_file, NamingConvention.UCSC)

        assert result.renamed_count == 3
        assert output_file.exists()

        content = output_file.read_text()
        assert ">chr1" in content
        assert ">chr2" in content
        assert ">chrX" in content


class TestCheckNcbiCompliance:
    """Tests for NCBI compliance checking."""

    def test_compliant_names(self):
        """Test compliant scaffold names."""
        scaffolds = [
            ("chr1", 1000, "ACGT"),
            ("scaffold_1", 2000, "GCTA"),
            ("contig.1", 3000, "TTTT"),
        ]
        result = check_ncbi_compliance(scaffolds)

        assert result.is_compliant

    def test_forbidden_characters(self):
        """Test detection of forbidden characters."""
        scaffolds = [
            ("chr 1", 1000, "ACGT"),  # Space
            ("chr>2", 2000, "GCTA"),  # >
            ("chr[3]", 3000, "TTTT"),  # Brackets
        ]
        result = check_ncbi_compliance(scaffolds)

        assert not result.is_compliant
        forbidden_issues = [i for i in result.issues if i.issue_type == "forbidden_characters"]
        assert len(forbidden_issues) == 3

    def test_duplicate_names(self):
        """Test detection of duplicate names."""
        scaffolds = [
            ("chr1", 1000, "ACGT"),
            ("chr2", 2000, "GCTA"),
            ("chr1", 3000, "TTTT"),  # Duplicate
        ]
        result = check_ncbi_compliance(scaffolds)

        assert not result.is_compliant
        dup_issues = [i for i in result.issues if i.issue_type == "duplicate_name"]
        assert len(dup_issues) == 1

    def test_long_names(self):
        """Test warning for long names."""
        long_name = "a" * 60
        scaffolds = [(long_name, 1000, "ACGT")]
        result = check_ncbi_compliance(scaffolds)

        # Long names are warnings, not errors
        assert result.is_compliant
        long_issues = [i for i in result.issues if i.issue_type == "long_name"]
        assert len(long_issues) == 1
        assert long_issues[0].severity == "warning"

    def test_non_standard_characters(self):
        """Test warning for non-standard characters."""
        scaffolds = [("chr1@special", 1000, "ACGT")]
        result = check_ncbi_compliance(scaffolds)

        # Non-standard chars are warnings
        assert result.is_compliant
        char_issues = [i for i in result.issues if i.issue_type == "non_standard_characters"]
        assert len(char_issues) == 1

    def test_just_names(self):
        """Test compliance check with just names (not tuples)."""
        names = ["chr1", "chr2", "chr 3"]
        result = check_ncbi_compliance(names)

        assert not result.is_compliant


class TestFormatFunctions:
    """Tests for formatting functions."""

    def test_format_compliance_summary_compliant(self):
        """Test formatting compliant result."""
        result = ComplianceResult(is_compliant=True, total_scaffolds=10)
        summary = format_compliance_summary(result)

        assert "NCBI SUBMISSION COMPLIANCE" in summary
        assert "COMPLIANT" in summary
        assert "10" in summary

    def test_format_compliance_summary_not_compliant(self):
        """Test formatting non-compliant result."""
        from chromdetect.standardize import ComplianceIssue
        result = ComplianceResult(
            is_compliant=False,
            total_scaffolds=10,
            issues=[
                ComplianceIssue("error", "forbidden_characters", "chr 1", "Has space"),
            ],
        )
        summary = format_compliance_summary(result)

        assert "NOT COMPLIANT" in summary
        assert "ERRORS:" in summary

    def test_format_rename_summary(self):
        """Test formatting rename result."""
        result = RenameResult(
            original_names=["chr1", "chr2"],
            new_names=["1", "2"],
            mapping={"chr1": "1", "chr2": "2"},
            unchanged=[],
            renamed_count=2,
        )
        summary = format_rename_summary(result)

        assert "SCAFFOLD RENAMING" in summary
        assert "Renamed:          2" in summary
        assert "chr1 -> 1" in summary


class TestRenameResultToDict:
    """Tests for RenameResult serialization."""

    def test_to_dict(self):
        """Test serialization to dictionary."""
        result = RenameResult(
            original_names=["chr1"],
            new_names=["1"],
            mapping={"chr1": "1"},
            unchanged=[],
            renamed_count=1,
        )
        d = result.to_dict()

        assert d["renamed_count"] == 1
        assert d["total_scaffolds"] == 1
        assert "chr1" in d["mapping"]


class TestStandardizationCLI:
    """Tests for standardization CLI commands."""

    @pytest.fixture
    def ucsc_fasta(self) -> Path:
        """Create a FASTA with UCSC naming."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">chr1\nACGTACGT\n")
            f.write(">chr2\nGCTAGCTA\n")
            f.write(">chrX\nTTTTAAAA\n")
            f.write(">chrM\nCCCCGGGG\n")
            f.flush()
            return Path(f.name)

    @pytest.fixture
    def ensembl_fasta(self) -> Path:
        """Create a FASTA with Ensembl naming."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">1\nACGTACGT\n")
            f.write(">2\nGCTAGCTA\n")
            f.write(">X\nTTTTAAAA\n")
            f.write(">MT\nCCCCGGGG\n")
            f.flush()
            return Path(f.name)

    @pytest.fixture
    def bad_fasta(self) -> Path:
        """Create a FASTA with NCBI-non-compliant names."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">chr 1\nACGTACGT\n")  # Space
            f.write(">chr>2\nGCTAGCTA\n")  # >
            f.write(">chr1\nTTTTAAAA\n")
            f.write(">chr1\nCCCCGGGG\n")  # Duplicate
            f.flush()
            return Path(f.name)

    def test_detect_convention(self, ucsc_fasta):
        """Test --detect-convention flag."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(ucsc_fasta),
                "--detect-convention",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "ucsc" in result.stdout.lower()

    def test_detect_convention_json(self, ensembl_fasta):
        """Test --detect-convention with JSON output."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(ensembl_fasta),
                "--detect-convention",
                "-f", "json",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert data["detected_convention"] == "ensembl"

    def test_check_compliance_pass(self, ucsc_fasta):
        """Test --check-compliance with compliant file."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(ucsc_fasta),
                "--check-compliance",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "COMPLIANT" in result.stdout

    def test_check_compliance_fail(self, bad_fasta):
        """Test --check-compliance with non-compliant file."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(bad_fasta),
                "--check-compliance",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "NOT COMPLIANT" in result.stdout

    def test_check_compliance_json(self, bad_fasta):
        """Test --check-compliance with JSON output."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(bad_fasta),
                "--check-compliance",
                "-f", "json",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        data = json.loads(result.stdout)
        assert data["is_compliant"] is False
        assert data["error_count"] > 0

    def test_rename_to_ensembl(self, ucsc_fasta):
        """Test --rename to ensembl convention."""
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as out:
            result = subprocess.run(
                [
                    sys.executable, "-m", "chromdetect",
                    str(ucsc_fasta),
                    "--rename", "ensembl",
                    "-o", out.name,
                    "-q",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0

            content = Path(out.name).read_text()
            assert ">1\n" in content
            assert ">2\n" in content
            assert ">X\n" in content
            assert ">MT\n" in content

    def test_rename_to_ucsc(self, ensembl_fasta):
        """Test --rename to UCSC convention."""
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as out:
            result = subprocess.run(
                [
                    sys.executable, "-m", "chromdetect",
                    str(ensembl_fasta),
                    "--rename", "ucsc",
                    "-o", out.name,
                    "-q",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0

            content = Path(out.name).read_text()
            assert ">chr1\n" in content
            assert ">chr2\n" in content
            assert ">chrX\n" in content
            assert ">chrM\n" in content

    def test_rename_to_refseq(self, ensembl_fasta):
        """Test --rename to RefSeq convention."""
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as out:
            result = subprocess.run(
                [
                    sys.executable, "-m", "chromdetect",
                    str(ensembl_fasta),
                    "--rename", "refseq",
                    "--species", "human",
                    "-o", out.name,
                    "-q",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0

            content = Path(out.name).read_text()
            assert "NC_000001" in content
            assert "NC_000002" in content

    def test_rename_json_output(self, ucsc_fasta):
        """Test --rename with JSON output."""
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as out:
            result = subprocess.run(
                [
                    sys.executable, "-m", "chromdetect",
                    str(ucsc_fasta),
                    "--rename", "ensembl",
                    "-o", out.name,
                    "-f", "json",
                    "-q",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0
            data = json.loads(result.stdout)
            assert data["renamed_count"] == 4

    def test_rename_requires_output(self, ucsc_fasta):
        """Test that --rename requires --output."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(ucsc_fasta),
                "--rename", "ensembl",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "requires --output" in result.stderr

    def test_help_includes_standardization(self):
        """Test that help text includes standardization options."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "--help"],
            capture_output=True,
            text=True,
        )
        assert "--rename" in result.stdout
        assert "--check-compliance" in result.stdout
        assert "--detect-convention" in result.stdout
