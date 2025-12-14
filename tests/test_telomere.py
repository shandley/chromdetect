"""
Tests for telomere detection functionality.
"""


from chromdetect.telomere import (
    TELOMERE_MOTIFS,
    TelomereResult,
    _count_repeats,
    detect_telomere,
    detect_telomeres_batch,
    get_telomere_summary,
)


class TestCountRepeats:
    """Tests for the _count_repeats helper function."""

    def test_count_vertebrate_motif(self):
        """Test counting TTAGGG repeats."""
        seq = "TTAGGGTTAGGGTTAGGG"  # 3 repeats
        assert _count_repeats(seq, "TTAGGG") == 3

    def test_count_no_repeats(self):
        """Test with no repeats."""
        seq = "ACGTACGTACGT"
        assert _count_repeats(seq, "TTAGGG") == 0

    def test_count_single_repeat(self):
        """Test with single repeat."""
        seq = "TTAGGG"
        assert _count_repeats(seq, "TTAGGG") == 1

    def test_count_case_insensitive(self):
        """Test case insensitivity."""
        seq = "ttagggttagggttaggg"
        assert _count_repeats(seq, "TTAGGG") == 3

    def test_count_empty_sequence(self):
        """Test with empty sequence."""
        assert _count_repeats("", "TTAGGG") == 0

    def test_count_empty_motif(self):
        """Test with empty motif."""
        assert _count_repeats("ACGT", "") == 0

    def test_count_multiple_regions(self):
        """Test with multiple separate repeat regions."""
        # Two regions of repeats separated by other sequence
        seq = "TTAGGGTTAGGGTTAGGG" + "ACGT" * 10 + "TTAGGGTTAGGG"
        # Should return the max (3)
        assert _count_repeats(seq, "TTAGGG") == 3


class TestDetectTelomere:
    """Tests for the detect_telomere function."""

    def test_detect_vertebrate_5prime(self):
        """Test detection at 5' end (reverse complement)."""
        # CCCTAA repeats at 5' end
        seq = "CCCTAACCCTAACCCTAACCCTAA" + "A" * 1000
        result = detect_telomere(seq)
        assert result.has_5prime
        assert not result.has_3prime
        assert result.motif_5prime == "CCCTAA"
        assert result.repeats_5prime >= 3

    def test_detect_vertebrate_3prime(self):
        """Test detection at 3' end (forward motif)."""
        # TTAGGG repeats at 3' end
        seq = "A" * 1000 + "TTAGGGTTAGGGTTAGGGTTAGGG"
        result = detect_telomere(seq)
        assert not result.has_5prime
        assert result.has_3prime
        assert result.motif_3prime == "TTAGGG"
        assert result.repeats_3prime >= 3

    def test_detect_both_ends(self):
        """Test T2T detection (both ends)."""
        seq = "CCCTAACCCTAACCCTAACCCTAA" + "A" * 1000 + "TTAGGGTTAGGGTTAGGGTTAGGG"
        result = detect_telomere(seq)
        assert result.has_5prime
        assert result.has_3prime
        assert result.is_complete

    def test_detect_no_telomeres(self):
        """Test sequence without telomeres."""
        seq = "ACGT" * 1000
        result = detect_telomere(seq)
        assert not result.has_5prime
        assert not result.has_3prime
        assert not result.has_telomere
        assert not result.is_complete

    def test_detect_empty_sequence(self):
        """Test with empty sequence."""
        result = detect_telomere("")
        assert not result.has_telomere

    def test_detect_plant_motif(self):
        """Test plant telomere motif (TTTAGGG)."""
        seq = "CCCTAAACCCTAAACCCTAAACCCTAAA" + "A" * 1000 + "TTTAGGGTTTAGGGTTTAGGGTTTAGGG"
        result = detect_telomere(seq)
        assert result.has_telomere
        # The strongest match should be detected
        assert result.organism_type is not None

    def test_detect_with_custom_window(self):
        """Test with custom search window."""
        # Create a sequence with telomere at the 3' end
        # Use a long sequence so the window matters
        telomere = "TTAGGG" * 10  # 60 bp telomere
        seq = "A" * 50000 + telomere  # 50,060 bp total

        # With 5000 bp window, we should find the telomere (it's in last 60 bp)
        result = detect_telomere(seq, search_window=5000)
        assert result.has_3prime

        # With 10 bp window, we won't find enough repeats
        result = detect_telomere(seq, search_window=10)
        # Window will be min(10, 50060//2) = 10, which only captures part of telomere
        assert result.repeats_3prime < 3  # Not enough repeats to count

    def test_detect_min_repeats(self):
        """Test minimum repeats threshold."""
        # Only 2 repeats - below default threshold
        seq = "A" * 1000 + "TTAGGGTTAGGG"
        result = detect_telomere(seq, min_repeats=3)
        assert not result.has_3prime

        # With lower threshold, should be detected
        result = detect_telomere(seq, min_repeats=2)
        assert result.has_3prime


class TestTelomereResult:
    """Tests for TelomereResult dataclass."""

    def test_has_telomere_property(self):
        """Test has_telomere property."""
        result = TelomereResult(has_5prime=True)
        assert result.has_telomere

        result = TelomereResult(has_3prime=True)
        assert result.has_telomere

        result = TelomereResult()
        assert not result.has_telomere

    def test_is_complete_property(self):
        """Test is_complete (T2T) property."""
        result = TelomereResult(has_5prime=True, has_3prime=True)
        assert result.is_complete

        result = TelomereResult(has_5prime=True)
        assert not result.is_complete

    def test_to_dict(self):
        """Test dictionary conversion."""
        result = TelomereResult(
            has_5prime=True,
            has_3prime=True,
            motif_5prime="CCCTAA",
            motif_3prime="TTAGGG",
            repeats_5prime=10,
            repeats_3prime=8,
            organism_type="vertebrate",
        )
        d = result.to_dict()
        assert d["has_5prime"] is True
        assert d["has_3prime"] is True
        assert d["has_telomere"] is True
        assert d["is_complete"] is True
        assert d["motif_5prime"] == "CCCTAA"
        assert d["organism_type"] == "vertebrate"


class TestDetectTelomeresBatch:
    """Tests for batch telomere detection."""

    def test_batch_detection(self):
        """Test batch detection on multiple scaffolds."""
        scaffolds = [
            ("chr1", 2000, "CCCTAACCCTAACCCTAACCCTAA" + "A" * 1000 + "TTAGGGTTAGGGTTAGGGTTAGGG"),
            ("chr2", 1000, "A" * 1000),  # No telomeres
            ("chr3", 1500, "A" * 1000 + "TTAGGGTTAGGGTTAGGGTTAGGG"),  # 3' only
        ]
        results = detect_telomeres_batch(scaffolds)

        assert len(results) == 3
        assert results["chr1"].is_complete
        assert not results["chr2"].has_telomere
        assert results["chr3"].has_3prime and not results["chr3"].has_5prime


class TestGetTelomereSummary:
    """Tests for telomere summary statistics."""

    def test_summary_statistics(self):
        """Test summary calculation."""
        results = {
            "chr1": TelomereResult(has_5prime=True, has_3prime=True),
            "chr2": TelomereResult(has_5prime=True, has_3prime=False),
            "chr3": TelomereResult(has_5prime=False, has_3prime=True),
            "chr4": TelomereResult(has_5prime=False, has_3prime=False),
        }
        summary = get_telomere_summary(results)

        assert summary["total_scaffolds"] == 4
        assert summary["with_telomere"] == 3
        assert summary["with_5prime"] == 2
        assert summary["with_3prime"] == 2
        assert summary["complete_t2t"] == 1
        assert summary["percent_with_telomere"] == 75.0
        assert summary["percent_complete"] == 25.0

    def test_summary_empty(self):
        """Test summary with no results."""
        summary = get_telomere_summary({})
        assert summary["total_scaffolds"] == 0
        assert summary["percent_with_telomere"] == 0.0


class TestTelomereMotifs:
    """Tests for telomere motif definitions."""

    def test_motifs_defined(self):
        """Test that telomere motifs are properly defined."""
        assert len(TELOMERE_MOTIFS) > 0

        for forward, reverse, name in TELOMERE_MOTIFS:
            # Forward and reverse should be proper DNA sequences
            assert all(c in "ACGT" for c in forward)
            assert all(c in "ACGT" for c in reverse)
            # Name should be non-empty
            assert name
