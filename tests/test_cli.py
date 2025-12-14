"""Tests for command-line interface."""

import json
import tempfile
from pathlib import Path

import pytest

from chromdetect.cli import format_output, main
from chromdetect.core import ScaffoldInfo, AssemblyStats, classify_scaffolds


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


class TestCLIIntegration:
    """Integration tests for CLI."""

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

    def test_classify_sample_fasta(self, sample_fasta: Path) -> None:
        """Test classifying sample FASTA file."""
        from chromdetect.core import parse_fasta, classify_scaffolds

        scaffolds = parse_fasta(sample_fasta)
        results, stats = classify_scaffolds(scaffolds)

        assert stats.total_scaffolds == 3
        assert stats.chromosome_count == 2
        assert stats.unplaced_count == 1
