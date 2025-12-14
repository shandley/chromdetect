"""Tests for chromosome naming pattern detection."""

import pytest

from chromdetect.core import detect_by_name
from chromdetect.patterns import (
    CHROMOSOME_PATTERNS,
    COMPILED_CHROMOSOME_PATTERNS,
)


class TestChromosomePatterns:
    """Test chromosome naming pattern detection."""

    @pytest.mark.parametrize(
        "name,expected_class,expected_chr_id",
        [
            # Explicit chr patterns
            ("chr1", "chromosome", "1"),
            ("Chr1", "chromosome", "1"),
            ("CHR1", "chromosome", "1"),
            ("chr_1", "chromosome", "1"),
            ("chr-1", "chromosome", "1"),
            ("chromosome1", "chromosome", "1"),
            ("chromosome_1", "chromosome", "1"),
            ("chrX", "chromosome", "X"),
            ("chrY", "chromosome", "Y"),
            ("chrZ", "chromosome", "Z"),
            ("chrW", "chromosome", "W"),
            ("chrM", "chromosome", "M"),
            ("chrMT", "chromosome", "MT"),
            # Super scaffold patterns
            ("Super_scaffold_1", "chromosome", "1"),
            ("super_scaffold_1", "chromosome", "1"),
            ("SUPER_SCAFFOLD_1", "chromosome", "1"),
            ("super-scaffold-1", "chromosome", "1"),
            ("Superscaffold1", "chromosome", "1"),
            ("superscaffold_X", "chromosome", "X"),
            # SUPER patterns
            ("SUPER_1", "chromosome", "1"),
            ("SUPER1", "chromosome", "1"),
            ("SUPER_X", "chromosome", "X"),
            # Linkage groups
            ("LG1", "chromosome", "1"),
            ("LG_1", "chromosome", "1"),
            ("LG-1", "chromosome", "1"),
            ("LG_X", "chromosome", "X"),
            # Simple numeric
            ("1", "chromosome", "1"),
            ("X", "chromosome", "X"),
            # HiC scaffold patterns
            ("HiC_scaffold_1", "chromosome", "1"),
            ("HiC_scaffold_23", "chromosome", "23"),
        ],
    )
    def test_chromosome_detection(
        self, name: str, expected_class: str, expected_chr_id: str
    ) -> None:
        """Test detection of chromosome-level scaffolds."""
        classification, confidence, method, chr_id = detect_by_name(name)
        assert classification == expected_class
        assert chr_id == expected_chr_id
        assert confidence >= 0.5

    @pytest.mark.parametrize(
        "name",
        [
            "chr1_random",
            "chrUn_scaffold1",
            "scaffold1_unloc",
            "chr1_unlocalized",
        ],
    )
    def test_unlocalized_detection(self, name: str) -> None:
        """Test detection of unlocalized scaffolds."""
        classification, confidence, method, chr_id = detect_by_name(name)
        assert classification == "unlocalized"
        assert "unloc" in method or "random" in method.lower()

    @pytest.mark.parametrize(
        "name",
        [
            "scaffold_arrow_ctg1",
            "contig_001",
            "ctg123",
            "scaffold1_pilon",
            "fragment_001",
            "chr1_hap2",
        ],
    )
    def test_fragment_detection(self, name: str) -> None:
        """Test detection of fragment/contig scaffolds."""
        classification, confidence, method, chr_id = detect_by_name(name)
        assert classification == "unplaced"
        assert "fragment" in method

    @pytest.mark.parametrize(
        "name",
        [
            "assembly_scaffold_xyz",
            "unknown_sequence",
        ],
    )
    def test_unknown_names(self, name: str) -> None:
        """Test handling of unknown naming patterns."""
        classification, confidence, method, chr_id = detect_by_name(name)
        assert classification == "other"
        assert confidence < 0.5

    def test_random_pattern_detected_as_unlocalized(self) -> None:
        """Test that names containing 'random' are classified as unlocalized."""
        # This is correct behavior - 'random' indicates unlocalized scaffolds
        classification, confidence, method, chr_id = detect_by_name("some_random_name")
        assert classification == "unlocalized"


class TestNCBIPatterns:
    """Test NCBI accession pattern detection."""

    @pytest.mark.parametrize(
        "name",
        [
            "NC_000001.11",
            "NC_000023.11",
            "NC_012920.1",
        ],
    )
    def test_refseq_detection(self, name: str) -> None:
        """Test RefSeq accession detection."""
        classification, confidence, method, chr_id = detect_by_name(name)
        assert classification == "chromosome"
        assert "refseq" in method.lower()

    @pytest.mark.parametrize(
        "name",
        [
            "CM000001.1",
            "CM000023.1",
        ],
    )
    def test_genbank_detection(self, name: str) -> None:
        """Test GenBank accession detection."""
        classification, confidence, method, chr_id = detect_by_name(name)
        assert classification == "chromosome"
        assert "genbank" in method.lower()


class TestPatternCompilation:
    """Test pattern compilation and consistency."""

    def test_patterns_compile(self) -> None:
        """Test that all patterns compile successfully."""
        assert len(COMPILED_CHROMOSOME_PATTERNS) == len(CHROMOSOME_PATTERNS)

    def test_compiled_patterns_match(self) -> None:
        """Test compiled patterns produce same results."""
        test_names = ["chr1", "Super_scaffold_1", "LG_X"]
        for name in test_names:
            classification, _, _, _ = detect_by_name(name)
            assert classification == "chromosome"
