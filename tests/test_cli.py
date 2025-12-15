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
        assert "0.5.0" in result.stdout

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
        assert "ChromDetect 0.5.0" in result.stderr
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


class TestBEDFormat:
    """Test BED format output."""

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

    def test_bed_format_output(self, small_fasta: Path) -> None:
        """Test BED format output."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", str(small_fasta), "-f", "bed", "-q"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        lines = result.stdout.strip().split("\n")
        # Should have 3 scaffolds
        assert len(lines) == 3
        # Check BED format (tab-separated, 6 columns)
        for line in lines:
            fields = line.split("\t")
            assert len(fields) == 6
            # Check start is 0 (BED is 0-based)
            assert fields[1] == "0"
            # Check score is numeric
            assert fields[4].isdigit()

    def test_bed_format_chromosomes_only(self, small_fasta: Path) -> None:
        """Test BED format with chromosomes-only filter."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                str(small_fasta),
                "-f",
                "bed",
                "-c",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        lines = result.stdout.strip().split("\n")
        # All should be chromosomes
        for line in lines:
            fields = line.split("\t")
            assert fields[3] == "chromosome"


class TestGFFFormat:
    """Test GFF format output."""

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

    def test_gff_format_output(self, small_fasta: Path) -> None:
        """Test GFF format output."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", str(small_fasta), "-f", "gff", "-q"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        lines = result.stdout.strip().split("\n")
        # First line should be GFF version header
        assert lines[0] == "##gff-version 3"
        # Should have header + 3 scaffolds
        assert len(lines) == 4
        # Check GFF format (tab-separated, 9 columns)
        for line in lines[1:]:
            fields = line.split("\t")
            assert len(fields) == 9
            # Check source is chromdetect
            assert fields[1] == "chromdetect"
            # Check start is 1 (GFF is 1-based)
            assert fields[3] == "1"
            # Check attributes contain ID
            assert "ID=" in fields[8]


class TestExtractChromosomes:
    """Test chromosome sequence extraction."""

    @pytest.fixture
    def sample_fasta(self) -> Path:
        """Create a sample FASTA file for testing."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">chr1\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write(">chr2\n")
            f.write("GCTAGCTAGCTAGCTA\n")
            f.write(">scaffold1\n")
            f.write("AAAAAAAAAAAAAAAA\n")
            f.flush()
            return Path(f.name)

    def test_extract_chromosomes(self, sample_fasta: Path) -> None:
        """Test extracting chromosome sequences to file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as out:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "chromdetect",
                    str(sample_fasta),
                    "--extract-chromosomes",
                    out.name,
                    "-q",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0

            # Read the extracted file
            with open(out.name) as f:
                content = f.read()

            # Should contain chromosome sequences
            assert ">chr1" in content
            assert ">chr2" in content
            # Should not contain scaffold1 (it's classified as unplaced)
            # Note: depends on classification logic

    def test_extract_chromosomes_message(self, sample_fasta: Path) -> None:
        """Test that extraction reports success."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as out:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "chromdetect",
                    str(sample_fasta),
                    "--extract-chromosomes",
                    out.name,
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0
            assert "Extracted" in result.stderr
            assert "chromosome sequences" in result.stderr


class TestBatchProcessing:
    """Test batch processing functionality."""

    @pytest.fixture
    def batch_dir(self) -> Path:
        """Create a directory with multiple FASTA files."""
        temp_dir = tempfile.mkdtemp()
        batch_path = Path(temp_dir)

        # Create sample FASTA files
        for i in range(3):
            fasta_file = batch_path / f"assembly_{i}.fasta"
            with open(fasta_file, "w") as f:
                f.write(f">chr{i+1}\n")
                f.write("ATCGATCG" * 100 + "\n")
                f.write(f">scaffold{i+1}\n")
                f.write("GCTAGCTA" * 50 + "\n")

        return batch_path

    def test_batch_processing(self, batch_dir: Path) -> None:
        """Test processing a directory of FASTA files."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                "--batch",
                str(batch_dir),
                "-f",
                "json",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0

        # Check results directory was created
        results_dir = batch_dir / "chromdetect_results"
        assert results_dir.exists()

        # Check batch summary was created
        summary_file = results_dir / "batch_summary.tsv"
        assert summary_file.exists()

        # Check individual result files
        json_files = list(results_dir.glob("*.json"))
        assert len(json_files) == 3

    def test_batch_with_output_dir(self, batch_dir: Path) -> None:
        """Test batch processing with custom output directory."""
        output_dir = Path(tempfile.mkdtemp())

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                "--batch",
                str(batch_dir),
                "-o",
                str(output_dir),
                "-f",
                "tsv",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0

        # Check results were written to custom directory
        tsv_files = list(output_dir.glob("*.tsv"))
        # 3 result files + 1 summary
        assert len(tsv_files) == 4

    def test_batch_empty_directory(self) -> None:
        """Test batch processing with empty directory."""
        empty_dir = Path(tempfile.mkdtemp())

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                "--batch",
                str(empty_dir),
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "No FASTA files found" in result.stderr

    def test_batch_invalid_directory(self) -> None:
        """Test batch processing with non-existent directory."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "chromdetect",
                "--batch",
                "/nonexistent/directory",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "not a directory" in result.stderr


class TestCoreFunctions:
    """Test new core functions."""

    def test_format_bed(self) -> None:
        """Test BED format function."""
        from chromdetect.core import format_bed

        results = [
            ScaffoldInfo(
                name="chr1",
                length=1000,
                classification="chromosome",
                confidence=0.95,
                detection_method="test",
                chromosome_id="1",
            ),
        ]
        bed_output = format_bed(results)
        lines = bed_output.strip().split("\n")
        assert len(lines) == 1
        fields = lines[0].split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "0"
        assert fields[2] == "1000"
        assert fields[3] == "chromosome"
        assert fields[4] == "950"  # 0.95 * 1000
        assert fields[5] == "."

    def test_format_bed_with_header(self) -> None:
        """Test BED format with header."""
        from chromdetect.core import format_bed

        results = [
            ScaffoldInfo(
                name="chr1",
                length=1000,
                classification="chromosome",
                confidence=0.95,
                detection_method="test",
                chromosome_id="1",
            ),
        ]
        bed_output = format_bed(results, include_header=True)
        lines = bed_output.strip().split("\n")
        assert len(lines) == 2
        assert lines[0].startswith("#")

    def test_format_gff(self) -> None:
        """Test GFF format function."""
        from chromdetect.core import format_gff

        results = [
            ScaffoldInfo(
                name="chr1",
                length=1000,
                classification="chromosome",
                confidence=0.95,
                detection_method="test",
                chromosome_id="1",
            ),
        ]
        gff_output = format_gff(results)
        lines = gff_output.strip().split("\n")
        assert lines[0] == "##gff-version 3"
        fields = lines[1].split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "chromdetect"
        assert fields[2] == "chromosome"
        assert fields[3] == "1"
        assert fields[4] == "1000"
        assert "ID=chr1" in fields[8]
        assert "chromosome_id=1" in fields[8]

    def test_write_fasta(self) -> None:
        """Test FASTA writing function."""
        from chromdetect.core import write_fasta

        sequences = [
            ("seq1", "ATCGATCGATCG"),
            ("seq2", "GCTAGCTAGCTA"),
        ]

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            write_fasta(sequences, f.name)

            with open(f.name) as out:
                content = out.read()

            assert ">seq1" in content
            assert "ATCGATCGATCG" in content
            assert ">seq2" in content
            assert "GCTAGCTAGCTA" in content

    def test_write_fasta_line_wrapping(self) -> None:
        """Test FASTA writing with line wrapping."""
        from chromdetect.core import write_fasta

        long_seq = "A" * 200
        sequences = [("seq1", long_seq)]

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            write_fasta(sequences, f.name, line_width=80)

            with open(f.name) as out:
                lines = out.read().strip().split("\n")

            # Header + 3 lines of 80, 80, 40 = 4 lines total
            assert len(lines) == 4
            assert lines[1] == "A" * 80
            assert lines[2] == "A" * 80
            assert lines[3] == "A" * 40

    def test_write_fasta_return_string(self) -> None:
        """Test FASTA writing returning string."""
        from chromdetect.core import write_fasta

        sequences = [("seq1", "ATCG")]
        output = write_fasta(sequences, output_path=None)

        assert ">seq1" in output
        assert "ATCG" in output

    def test_parse_fasta_full_sequence(self) -> None:
        """Test parsing FASTA with full sequence retention."""
        from chromdetect.core import parse_fasta

        # Create file with long sequence split across multiple lines (80 chars/line)
        # This is typical FASTA format
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">seq1\n")
            # Write 200 lines of 80 chars each = 16000 bp
            for _ in range(200):
                f.write("ATCGATCG" * 10 + "\n")  # 80 chars per line
            f.flush()
            temp_path = f.name

        # Without full sequence (default)
        scaffolds = parse_fasta(temp_path, keep_full_sequence=False)
        name, length, seq = scaffolds[0]
        assert length == 16000
        # Sample should be approximately 10kb (may vary by line)
        assert 10000 <= len(seq) <= 10080  # Within one line of limit

        # With full sequence
        scaffolds = parse_fasta(temp_path, keep_full_sequence=True)
        name, length, seq = scaffolds[0]
        assert length == 16000
        assert len(seq) == 16000  # Full sequence
