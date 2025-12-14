"""Tests for command-line interface."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from chromdetect.cli import (
    EXIT_DATAERR,
    EXIT_NOINPUT,
    format_output,
    show_patterns,
)
from chromdetect.core import AssemblyStats, ScaffoldInfo


class TestFormatOutput:
    """Test output formatting."""

    @pytest.fixture
    def sample_data(self) -> tuple[list[ScaffoldInfo], AssemblyStats]:
        """Create sample data for formatting tests."""
        results = [
            ScaffoldInfo(
                name="chr1",
                length=100_000_000,
                classification="chromosome",
                confidence=0.95,
                detection_method="name_chr_explicit",
                chromosome_id="1",
            ),
            ScaffoldInfo(
                name="scaffold1",
                length=500_000,
                classification="unplaced",
                confidence=0.6,
                detection_method="size_small",
                chromosome_id=None,
            ),
        ]
        stats = AssemblyStats(
            total_scaffolds=2,
            total_length=100_500_000,
            n50=100_000_000,
            n90=500_000,
            chromosome_count=1,
            chromosome_length=100_000_000,
            chromosome_n50=100_000_000,
            unlocalized_count=0,
            unplaced_count=1,
            largest_scaffold=100_000_000,
            gc_content=0.42,
        )
        return results, stats

    def test_json_format(
        self, sample_data: tuple[list[ScaffoldInfo], AssemblyStats]
    ) -> None:
        """Test JSON output format."""
        results, stats = sample_data
        output = format_output(results, stats, "json")

        data = json.loads(output)
        assert "summary" in data
        assert "scaffolds" in data
        assert data["summary"]["total_scaffolds"] == 2
        assert len(data["scaffolds"]) == 2

    def test_tsv_format(
        self, sample_data: tuple[list[ScaffoldInfo], AssemblyStats]
    ) -> None:
        """Test TSV output format."""
        results, stats = sample_data
        output = format_output(results, stats, "tsv")

        lines = output.strip().split("\n")
        assert len(lines) == 3  # Header + 2 scaffolds
        assert "name\tlength\tclassification" in lines[0]
        assert "chr1" in lines[1]

    def test_summary_format(
        self, sample_data: tuple[list[ScaffoldInfo], AssemblyStats]
    ) -> None:
        """Test summary output format."""
        results, stats = sample_data
        output = format_output(results, stats, "summary")

        assert "CHROMDETECT ASSEMBLY ANALYSIS" in output
        assert "Total scaffolds" in output
        assert "Chromosomes" in output
        assert "chr1" in output

    def test_invalid_format(
        self, sample_data: tuple[list[ScaffoldInfo], AssemblyStats]
    ) -> None:
        """Test invalid format raises error."""
        results, stats = sample_data
        with pytest.raises(ValueError, match="Unknown format"):
            format_output(results, stats, "invalid")


class TestShowPatterns:
    """Test pattern listing functionality."""

    def test_show_patterns_output(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Test that show_patterns outputs pattern information."""
        show_patterns()
        captured = capsys.readouterr()

        assert "ChromDetect Supported Naming Patterns" in captured.out
        assert "CHROMOSOME PATTERNS" in captured.out
        assert "UNLOCALIZED PATTERNS" in captured.out
        assert "FRAGMENT PATTERNS" in captured.out
        assert "chr_explicit" in captured.out


class TestCLIIntegration:
    """Integration tests for CLI using subprocess."""

    @pytest.fixture
    def sample_fasta(self) -> Path:
        """Create a sample FASTA file for testing."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">chr1\n")
            f.write("A" * 50_000_000 + "\n")
            f.write(">chr2\n")
            f.write("G" * 40_000_000 + "\n")
            f.write(">scaffold_ctg1\n")
            f.write("C" * 100_000 + "\n")
            f.flush()
            return Path(f.name)

    @pytest.fixture
    def small_fasta(self) -> Path:
        """Create a small FASTA file for quick tests."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">chr1\n")
            f.write("ATCGATCG" * 100 + "\n")
            f.write(">chr2\n")
            f.write("GCTAGCTA" * 100 + "\n")
            f.write(">scaffold1\n")
            f.write("AAAAAAAA" * 50 + "\n")
            f.flush()
            return Path(f.name)

    def test_classify_sample_fasta(self, sample_fasta: Path) -> None:
        """Test classifying sample FASTA file."""
        from chromdetect.core import classify_scaffolds, parse_fasta

        scaffolds = parse_fasta(sample_fasta)
        results, stats = classify_scaffolds(scaffolds)

        assert stats.total_scaffolds == 3
        assert stats.chromosome_count == 2
        assert stats.unplaced_count == 1

    def test_version_flag(self) -> None:
        """Test --version flag."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "--version"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "chromdetect" in result.stdout
        assert "0.2.0" in result.stdout

    def test_help_flag(self) -> None:
        """Test --help flag."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Detect chromosome-level scaffolds" in result.stdout
        assert "--format" in result.stdout
        assert "--karyotype" in result.stdout

    def test_list_patterns_flag(self) -> None:
        """Test --list-patterns flag."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "--list-patterns"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "CHROMOSOME PATTERNS" in result.stdout
        assert "chr_explicit" in result.stdout

    def test_json_format_output(self, small_fasta: Path) -> None:
        """Test JSON format output."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", str(small_fasta), "-f", "json", "-q"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert "summary" in data
        assert "scaffolds" in data
        assert data["summary"]["total_scaffolds"] == 3

    def test_tsv_format_output(self, small_fasta: Path) -> None:
        """Test TSV format output."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", str(small_fasta), "-f", "tsv", "-q"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        lines = result.stdout.strip().split("\n")
        assert len(lines) == 4  # Header + 3 scaffolds
        assert "name\tlength" in lines[0]

    def test_summary_format_output(self, small_fasta: Path) -> None:
        """Test summary format output."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                str(small_fasta),
                "-f",
                "summary",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "CHROMDETECT ASSEMBLY ANALYSIS" in result.stdout

    def test_quiet_flag(self, small_fasta: Path) -> None:
        """Test --quiet flag suppresses progress messages."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", str(small_fasta), "-q", "-f", "json"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Parsing" not in result.stderr
        assert "Found" not in result.stderr

    def test_verbose_flag(self, small_fasta: Path) -> None:
        """Test --verbose flag shows detailed info."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", str(small_fasta), "-v", "-f", "json"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "ChromDetect 0.2.0" in result.stderr
        assert "Input file:" in result.stderr or "Input:" in result.stderr

    def test_chromosomes_only_filter(self, small_fasta: Path) -> None:
        """Test --chromosomes-only filter."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                str(small_fasta),
                "-c",
                "-f",
                "json",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        # All remaining scaffolds should be chromosomes
        for scaffold in data["scaffolds"]:
            assert scaffold["classification"] == "chromosome"

    def test_min_confidence_filter(self, small_fasta: Path) -> None:
        """Test --min-confidence filter."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                str(small_fasta),
                "--min-confidence",
                "0.8",
                "-f",
                "json",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        for scaffold in data["scaffolds"]:
            assert scaffold["confidence"] >= 0.8

    def test_min_length_filter(self, small_fasta: Path) -> None:
        """Test --min-length filter."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                str(small_fasta),
                "--min-length",
                "500",
                "-f",
                "json",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        for scaffold in data["scaffolds"]:
            assert scaffold["length"] >= 500

    def test_output_file(self, small_fasta: Path) -> None:
        """Test --output flag writes to file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as out:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "chromdetect",
                    str(small_fasta),
                    "-f",
                    "json",
                    "-o",
                    out.name,
                    "-q",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0

            # Verify file was written
            with open(out.name) as f:
                data = json.load(f)
            assert "summary" in data

    def test_karyotype_option(self, small_fasta: Path) -> None:
        """Test --karyotype option."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                str(small_fasta),
                "-k",
                "2",
                "-f",
                "json",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        # Should adjust to match expected karyotype
        assert data["summary"]["chromosome_count"] <= 3


class TestCLIErrorHandling:
    """Test CLI error handling."""

    def test_missing_file_error(self) -> None:
        """Test error when file doesn't exist."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "/nonexistent/file.fasta"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == EXIT_NOINPUT
        assert "File not found" in result.stderr

    def test_invalid_fasta_error(self) -> None:
        """Test error for invalid FASTA format."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            f.write("This is not a FASTA file\n")
            f.write("Just some random text\n")
            f.flush()

            result = subprocess.run(
                [sys.executable, "-m", "chromdetect", f.name, "-q"],
                capture_output=True,
                text=True,
            )
            assert result.returncode == EXIT_DATAERR
            assert "Invalid FASTA format" in result.stderr

    def test_empty_file_error(self) -> None:
        """Test error for empty file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.flush()  # Empty file

            result = subprocess.run(
                [sys.executable, "-m", "chromdetect", f.name, "-q"],
                capture_output=True,
                text=True,
            )
            assert result.returncode == EXIT_DATAERR

    def test_missing_argument_error(self) -> None:
        """Test error when no file argument provided."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect"],
            capture_output=True,
            text=True,
        )
        # Should exit with usage error
        assert result.returncode != 0
        assert "required: fasta" in result.stderr or "required" in result.stderr


class TestStdinSupport:
    """Test stdin input support."""

    def test_stdin_input(self) -> None:
        """Test reading from stdin with '-'."""
        fasta_content = ">chr1\nATCGATCG\n>chr2\nGCTAGCTA\n"
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "-", "-f", "json", "-q"],
            input=fasta_content,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert data["summary"]["total_scaffolds"] == 2

    def test_stdin_empty_error(self) -> None:
        """Test error for empty stdin."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "-", "-q"],
            input="",
            capture_output=True,
            text=True,
        )
        assert result.returncode == EXIT_DATAERR

    def test_stdin_invalid_fasta_error(self) -> None:
        """Test error for invalid FASTA from stdin."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "-", "-q"],
            input="not a fasta file\n",
            capture_output=True,
            text=True,
        )
        assert result.returncode == EXIT_DATAERR
        assert "Invalid FASTA format" in result.stderr
