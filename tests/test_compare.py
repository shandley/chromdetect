"""
Tests for assembly comparison functionality.
"""


import pytest

from chromdetect.compare import (
    ComparisonResult,
    compare_assemblies,
    format_comparison_summary,
    format_comparison_tsv,
)
from chromdetect.core import AssemblyStats, ScaffoldInfo


def make_scaffold(name, length, classification, confidence=0.9, chr_id=None):
    """Helper to create ScaffoldInfo for testing."""
    return ScaffoldInfo(
        name=name,
        length=length,
        classification=classification,
        confidence=confidence,
        detection_method="test",
        chromosome_id=chr_id,
    )


def make_stats(
    total_scaffolds=10,
    total_length=1000000,
    n50=100000,
    n90=50000,
    chromosome_count=5,
    chromosome_length=800000,
    chromosome_n50=160000,
    unlocalized_count=2,
    unplaced_count=3,
    largest_scaffold=200000,
    gc_content=0.4,
):
    """Helper to create AssemblyStats for testing."""
    return AssemblyStats(
        total_scaffolds=total_scaffolds,
        total_length=total_length,
        n50=n50,
        n90=n90,
        chromosome_count=chromosome_count,
        chromosome_length=chromosome_length,
        chromosome_n50=chromosome_n50,
        unlocalized_count=unlocalized_count,
        unplaced_count=unplaced_count,
        largest_scaffold=largest_scaffold,
        gc_content=gc_content,
    )


class TestCompareAssemblies:
    """Tests for compare_assemblies function."""

    def test_identical_assemblies(self):
        """Test comparing identical assemblies."""
        scaffolds = [
            make_scaffold("chr1", 100000, "chromosome", chr_id="1"),
            make_scaffold("chr2", 80000, "chromosome", chr_id="2"),
            make_scaffold("scaffold1", 10000, "unplaced"),
        ]
        stats = make_stats()

        result = compare_assemblies(
            scaffolds, stats, scaffolds, stats,
            "Assembly1", "Assembly2"
        )

        assert len(result.shared_chromosomes) == 2
        assert len(result.unique_to_1) == 0
        assert len(result.unique_to_2) == 0
        assert len(result.size_differences) == 0

    def test_different_assemblies(self):
        """Test comparing different assemblies."""
        scaffolds1 = [
            make_scaffold("chr1", 100000, "chromosome", chr_id="1"),
            make_scaffold("chr2", 80000, "chromosome", chr_id="2"),
            make_scaffold("chrX", 90000, "chromosome", chr_id="X"),
        ]
        scaffolds2 = [
            make_scaffold("chr1", 105000, "chromosome", chr_id="1"),
            make_scaffold("chr2", 80000, "chromosome", chr_id="2"),
            make_scaffold("chr3", 70000, "chromosome", chr_id="3"),
        ]
        stats = make_stats()

        result = compare_assemblies(
            scaffolds1, stats, scaffolds2, stats,
            "Assembly1", "Assembly2"
        )

        assert "chr1" in result.shared_chromosomes
        assert "chr2" in result.shared_chromosomes
        assert "chrX" in result.unique_to_1
        assert "chr3" in result.unique_to_2
        # chr1 has different size
        assert "chr1" in result.size_differences
        assert result.size_differences["chr1"] == 5000  # 105000 - 100000

    def test_classification_changes(self):
        """Test detecting classification changes."""
        scaffolds1 = [
            make_scaffold("scaffold1", 100000, "unplaced"),
            make_scaffold("chr1", 80000, "chromosome"),
        ]
        scaffolds2 = [
            make_scaffold("scaffold1", 100000, "chromosome"),  # Changed!
            make_scaffold("chr1", 80000, "unlocalized"),  # Changed!
        ]
        stats = make_stats()

        result = compare_assemblies(
            scaffolds1, stats, scaffolds2, stats,
            "Assembly1", "Assembly2"
        )

        assert len(result.classification_changes) == 2
        changes = {name: (c1, c2) for name, c1, c2 in result.classification_changes}
        assert changes["scaffold1"] == ("unplaced", "chromosome")
        assert changes["chr1"] == ("chromosome", "unlocalized")

    def test_summary_calculation(self):
        """Test summary statistics calculation."""
        stats1 = make_stats(n50=100000, chromosome_count=5)
        stats2 = make_stats(n50=150000, chromosome_count=6)

        result = compare_assemblies(
            [], stats1, [], stats2,
            "Assembly1", "Assembly2"
        )

        summary = result.summary()
        assert summary["n50_difference"] == 50000
        assert summary["chromosome_count_difference"] == 1


class TestComparisonResultToDict:
    """Tests for ComparisonResult.to_dict method."""

    def test_to_dict_structure(self):
        """Test that to_dict produces expected structure."""
        stats = make_stats()
        result = ComparisonResult(
            assembly1_name="Asm1",
            assembly2_name="Asm2",
            stats1=stats,
            stats2=stats,
            shared_chromosomes=["chr1", "chr2"],
            unique_to_1=["chrX"],
            unique_to_2=["chr3"],
            size_differences={"chr1": 1000},
            classification_changes=[("scaffold1", "unplaced", "chromosome")],
        )

        d = result.to_dict()

        assert d["assembly1_name"] == "Asm1"
        assert d["assembly2_name"] == "Asm2"
        assert "stats1" in d
        assert "stats2" in d
        assert d["shared_chromosomes"] == ["chr1", "chr2"]
        assert d["unique_to_1"] == ["chrX"]
        assert d["unique_to_2"] == ["chr3"]
        assert "summary" in d


class TestFormatComparisonSummary:
    """Tests for format_comparison_summary function."""

    def test_summary_format(self):
        """Test that summary format contains expected sections."""
        stats = make_stats()
        result = ComparisonResult(
            assembly1_name="TestAsm1",
            assembly2_name="TestAsm2",
            stats1=stats,
            stats2=stats,
            shared_chromosomes=["chr1"],
            unique_to_1=[],
            unique_to_2=[],
            size_differences={},
            classification_changes=[],
        )

        output = format_comparison_summary(result)

        assert "CHROMDETECT ASSEMBLY COMPARISON" in output
        assert "TestAsm1" in output
        assert "TestAsm2" in output
        assert "STATISTICS COMPARISON" in output
        assert "CHROMOSOME COMPARISON" in output
        assert "Shared chromosomes" in output

    def test_summary_with_differences(self):
        """Test summary with actual differences."""
        stats1 = make_stats(n50=100000)
        stats2 = make_stats(n50=200000)
        result = ComparisonResult(
            assembly1_name="Asm1",
            assembly2_name="Asm2",
            stats1=stats1,
            stats2=stats2,
            shared_chromosomes=["chr1"],
            unique_to_1=["chrX"],
            unique_to_2=["chr3", "chr4"],
            size_differences={"chr1": 5000},
            classification_changes=[("scaffold1", "unplaced", "chromosome")],
        )

        output = format_comparison_summary(result)

        assert "chrX" in output or "Unique to Assembly 1" in output
        assert "SIZE DIFFERENCES" in output
        assert "CLASSIFICATION CHANGES" in output


class TestFormatComparisonTSV:
    """Tests for format_comparison_tsv function."""

    def test_tsv_format(self):
        """Test that TSV format is valid."""
        stats = make_stats()
        result = ComparisonResult(
            assembly1_name="Asm1",
            assembly2_name="Asm2",
            stats1=stats,
            stats2=stats,
            shared_chromosomes=[],
            unique_to_1=[],
            unique_to_2=[],
            size_differences={},
            classification_changes=[],
        )

        output = format_comparison_tsv(result)
        lines = output.strip().split("\n")

        # Check header
        assert lines[0] == "metric\tassembly1\tassembly2\tdifference"

        # Check that each line has 4 columns
        for line in lines[1:]:
            columns = line.split("\t")
            assert len(columns) == 4

    def test_tsv_metrics(self):
        """Test that expected metrics are in TSV output."""
        stats = make_stats()
        result = ComparisonResult(
            assembly1_name="Asm1",
            assembly2_name="Asm2",
            stats1=stats,
            stats2=stats,
            shared_chromosomes=[],
            unique_to_1=[],
            unique_to_2=[],
            size_differences={},
            classification_changes=[],
        )

        output = format_comparison_tsv(result)

        assert "n50\t" in output
        assert "chromosome_count\t" in output
        assert "gc_content\t" in output
