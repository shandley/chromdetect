"""
Tests for convenience functions (classify_fasta and compare_fasta_files).

These functions provide a simpler API for common use cases.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from chromdetect import (
    classify_fasta,
    compare_fasta_files,
)


@pytest.fixture
def simple_fasta() -> str:
    """Create a simple test FASTA file."""
    content = """>chr1 test chromosome 1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2 test chromosome 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>scaffold_1 unplaced
NNNNNNNNNN
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(content)
        return f.name


@pytest.fixture
def second_fasta() -> str:
    """Create a second test FASTA file for comparison."""
    content = """>chr1 test chromosome 1 - improved
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2 test chromosome 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>chr3 new chromosome
ATATATATATATATATATATATATATATATATATATATATATATATATAT
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(content)
        return f.name


class TestClassifyFasta:
    """Tests for classify_fasta convenience function."""

    def test_basic_classification(self, simple_fasta: str) -> None:
        """Test basic FASTA classification."""
        results, stats = classify_fasta(simple_fasta)

        assert stats.total_scaffolds == 3
        assert stats.chromosome_count == 2
        assert stats.unplaced_count == 1

    def test_returns_scaffold_info(self, simple_fasta: str) -> None:
        """Test that results contain proper ScaffoldInfo objects."""
        results, stats = classify_fasta(simple_fasta)

        assert len(results) == 3
        chr1 = next(r for r in results if r.name == "chr1")
        assert chr1.classification == "chromosome"
        assert chr1.chromosome_id == "1"
        assert chr1.gc_content is not None

    def test_with_expected_chromosomes(self, simple_fasta: str) -> None:
        """Test with karyotype adjustment."""
        results, stats = classify_fasta(simple_fasta, expected_chromosomes=2)

        # Should find exactly 2 chromosomes
        assert stats.chromosome_count == 2

    def test_with_min_chromosome_size(self, simple_fasta: str) -> None:
        """Test with custom minimum chromosome size."""
        # Set very high min size - nothing should be classified as chromosome by size
        results, stats = classify_fasta(simple_fasta, min_chromosome_size=1_000_000_000)

        # Should still find chromosomes by name pattern
        assert stats.chromosome_count == 2

    def test_accepts_path_object(self, simple_fasta: str) -> None:
        """Test that Path objects are accepted."""
        results, stats = classify_fasta(Path(simple_fasta))

        assert stats.total_scaffolds == 3

    def test_gc_content_calculated(self, simple_fasta: str) -> None:
        """Test that GC content is calculated for scaffolds."""
        results, stats = classify_fasta(simple_fasta)

        # Check overall GC content
        assert stats.gc_content is not None

        # Check per-scaffold GC content
        chr1 = next(r for r in results if r.name == "chr1")
        assert chr1.gc_content is not None
        assert 0.0 <= chr1.gc_content <= 1.0


class TestCompareFastaFiles:
    """Tests for compare_fasta_files convenience function."""

    def test_basic_comparison(self, simple_fasta: str, second_fasta: str) -> None:
        """Test basic FASTA file comparison."""
        result = compare_fasta_files(simple_fasta, second_fasta)

        assert result.assembly1_name is not None
        assert result.assembly2_name is not None
        assert result.stats1.total_scaffolds == 3
        assert result.stats2.total_scaffolds == 3

    def test_shared_chromosomes(self, simple_fasta: str, second_fasta: str) -> None:
        """Test detection of shared chromosomes."""
        result = compare_fasta_files(simple_fasta, second_fasta)

        # chr1 and chr2 are in both
        assert "chr1" in result.shared_chromosomes
        assert "chr2" in result.shared_chromosomes

    def test_unique_chromosomes(self, simple_fasta: str, second_fasta: str) -> None:
        """Test detection of unique chromosomes."""
        result = compare_fasta_files(simple_fasta, second_fasta)

        # chr3 is only in second assembly
        assert "chr3" in result.unique_to_2

    def test_size_differences(self, simple_fasta: str, second_fasta: str) -> None:
        """Test detection of size differences."""
        result = compare_fasta_files(simple_fasta, second_fasta)

        # chr1 is larger in second assembly
        if "chr1" in result.size_differences:
            assert result.size_differences["chr1"] > 0

    def test_summary_method(self, simple_fasta: str, second_fasta: str) -> None:
        """Test summary method returns expected keys."""
        result = compare_fasta_files(simple_fasta, second_fasta)
        summary = result.summary()

        assert "total_shared_chromosomes" in summary
        assert "unique_to_1" in summary
        assert "unique_to_2" in summary
        assert "n50_difference" in summary
        assert "chromosome_count_difference" in summary

    def test_to_dict_method(self, simple_fasta: str, second_fasta: str) -> None:
        """Test to_dict serialization."""
        result = compare_fasta_files(simple_fasta, second_fasta)
        d = result.to_dict()

        assert "assembly1_name" in d
        assert "assembly2_name" in d
        assert "stats1" in d
        assert "stats2" in d
        assert "shared_chromosomes" in d

    def test_extracts_assembly_names(
        self, simple_fasta: str, second_fasta: str
    ) -> None:
        """Test that assembly names are extracted from file paths."""
        result = compare_fasta_files(simple_fasta, second_fasta)

        # Names should be derived from file stems
        assert result.assembly1_name != ""
        assert result.assembly2_name != ""
        # Should not include .fasta extension
        assert ".fasta" not in result.assembly1_name
        assert ".fasta" not in result.assembly2_name


class TestIntegration:
    """Integration tests for convenience functions."""

    def test_classify_then_compare_workflow(
        self, simple_fasta: str, second_fasta: str
    ) -> None:
        """Test a typical workflow using both convenience functions."""
        # First classify both assemblies
        results1, stats1 = classify_fasta(simple_fasta)
        results2, stats2 = classify_fasta(second_fasta)

        # Then compare them using the file-based convenience function
        comparison = compare_fasta_files(simple_fasta, second_fasta)

        # Results should be consistent
        assert comparison.stats1.chromosome_count == stats1.chromosome_count
        assert comparison.stats2.chromosome_count == stats2.chromosome_count

    def test_import_from_package_root(self) -> None:
        """Test that convenience functions can be imported from package root."""
        from chromdetect import classify_fasta, compare_fasta_files

        assert callable(classify_fasta)
        assert callable(compare_fasta_files)
