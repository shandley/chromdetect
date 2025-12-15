"""
Tests for GC content calculation functionality.
"""


import pytest

from chromdetect.core import calculate_gc, classify_scaffolds


class TestCalculateGC:
    """Tests for the calculate_gc function."""

    def test_empty_sequence(self):
        """Test GC calculation with empty sequence."""
        assert calculate_gc("") is None

    def test_all_gc(self):
        """Test sequence with all GC bases."""
        assert calculate_gc("GCGCGCGC") == 1.0

    def test_all_at(self):
        """Test sequence with all AT bases."""
        assert calculate_gc("ATATATAT") == 0.0

    def test_half_gc(self):
        """Test sequence with 50% GC content."""
        result = calculate_gc("GCATAT")
        assert result == pytest.approx(1 / 3, rel=0.01)

    def test_mixed_case(self):
        """Test that mixed case is handled."""
        assert calculate_gc("gcGC") == calculate_gc("GCGC")

    def test_with_n_bases(self):
        """Test that N bases are not counted."""
        # Only count A, T, G, C
        result = calculate_gc("GCNNNN")
        # GC = 2, total = 2 (N's not counted)
        assert result == 1.0

    def test_realistic_sequence(self):
        """Test with a realistic sequence."""
        # 60% GC content
        seq = "GCGCGCATAT"  # 6 GC, 4 AT
        result = calculate_gc(seq)
        assert result == pytest.approx(0.6, rel=0.01)


class TestGCContentInScaffoldInfo:
    """Tests for GC content in ScaffoldInfo."""

    def test_gc_content_calculated(self):
        """Test that GC content is calculated for scaffolds."""
        # Create a scaffold with known GC content
        # Sequence is 60% GC
        seq = "GCGCGCATAT" * 1000  # 10kb
        scaffolds = [("test_scaffold", len(seq), seq)]

        results, stats = classify_scaffolds(scaffolds)

        # Check per-scaffold GC
        assert len(results) == 1
        assert results[0].gc_content is not None
        assert results[0].gc_content == pytest.approx(0.6, rel=0.01)

    def test_gc_content_in_stats(self):
        """Test that overall GC content is in stats."""
        seq = "GCGCATATAT" * 1000  # 10kb, 40% GC
        scaffolds = [("test_scaffold", len(seq), seq)]

        results, stats = classify_scaffolds(scaffolds)

        assert stats.gc_content is not None
        assert stats.gc_content == pytest.approx(0.4, rel=0.01)

    def test_gc_content_none_for_empty_sequence(self):
        """Test that GC content is None when sequence is empty."""
        scaffolds = [("test_scaffold", 1000, "")]  # Empty sequence

        results, stats = classify_scaffolds(scaffolds)

        assert len(results) == 1
        assert results[0].gc_content is None

    def test_gc_content_to_dict(self):
        """Test that GC content is included in to_dict output."""
        seq = "GCGCGCATAT" * 1000
        scaffolds = [("test_scaffold", len(seq), seq)]

        results, stats = classify_scaffolds(scaffolds)

        result_dict = results[0].to_dict()
        assert "gc_content" in result_dict
        assert result_dict["gc_content"] == pytest.approx(0.6, rel=0.01)


class TestGCContentEdgeCases:
    """Edge case tests for GC content."""

    def test_very_short_sequence(self):
        """Test GC content with very short sequence."""
        assert calculate_gc("GC") == 1.0
        assert calculate_gc("AT") == 0.0
        assert calculate_gc("GCAT") == 0.5

    def test_ambiguous_bases(self):
        """Test handling of ambiguous bases (N, R, Y, etc.)."""
        # Only ACGT should be counted
        result = calculate_gc("GCRYNNNN")  # Only G, C count
        assert result == 1.0  # Both GC

    def test_lowercase_input(self):
        """Test that lowercase input works."""
        assert calculate_gc("gcgcatat") == calculate_gc("GCGCATAT")
