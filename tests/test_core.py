"""Tests for core classification functionality."""

import tempfile
from pathlib import Path

import pytest

from chromdetect.core import (
    AssemblyStats,
    ScaffoldInfo,
    calculate_gc,
    calculate_n50,
    calculate_n90,
    classify_scaffolds,
    detect_by_size,
    parse_fasta,
)


class TestN50Calculations:
    """Test N50/N90 calculations."""

    def test_n50_simple(self) -> None:
        """Test N50 with simple list."""
        lengths = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10]
        # Total = 550, 50% = 275
        # Running: 100->190->270->340 (crosses at 70)
        assert calculate_n50(lengths) == 70

    def test_n50_equal(self) -> None:
        """Test N50 with equal lengths."""
        lengths = [100, 100, 100, 100]
        assert calculate_n50(lengths) == 100

    def test_n50_empty(self) -> None:
        """Test N50 with empty list."""
        assert calculate_n50([]) == 0

    def test_n50_single(self) -> None:
        """Test N50 with single scaffold."""
        assert calculate_n50([1000]) == 1000

    def test_n90_simple(self) -> None:
        """Test N90 calculation."""
        lengths = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10]
        # Total = 550, 90% = 495
        # Running: 100->190->270->340->400->450->490->520 (crosses at 30)
        n90 = calculate_n90(lengths)
        assert n90 > 0
        assert n90 == 30  # N90 is the length at which we reach 90% of total


class TestGCCalculation:
    """Test GC content calculation."""

    def test_gc_50_percent(self) -> None:
        """Test 50% GC content."""
        assert calculate_gc("ATGC") == 0.5

    def test_gc_100_percent(self) -> None:
        """Test 100% GC content."""
        assert calculate_gc("GCGC") == 1.0

    def test_gc_0_percent(self) -> None:
        """Test 0% GC content."""
        assert calculate_gc("ATAT") == 0.0

    def test_gc_empty(self) -> None:
        """Test empty sequence."""
        assert calculate_gc("") is None

    def test_gc_case_insensitive(self) -> None:
        """Test case insensitivity."""
        assert calculate_gc("atgc") == calculate_gc("ATGC")


class TestSizeDetection:
    """Test size-based detection."""

    def test_large_scaffold(self) -> None:
        """Test detection of large scaffolds."""
        classification, confidence, method = detect_by_size(
            length=50_000_000,
            n50=30_000_000,
            largest=100_000_000,
        )
        assert classification == "chromosome"
        assert confidence >= 0.7
        assert "large" in method

    def test_n50_scaffold(self) -> None:
        """Test detection of N50-relative scaffolds."""
        classification, confidence, method = detect_by_size(
            length=20_000_000,
            n50=30_000_000,
            largest=100_000_000,
            min_chromosome_size=50_000_000,
        )
        assert classification == "chromosome"
        assert "n50" in method

    def test_small_scaffold(self) -> None:
        """Test detection of small scaffolds."""
        classification, confidence, method = detect_by_size(
            length=100_000,
            n50=30_000_000,
            largest=100_000_000,
        )
        assert classification == "unplaced"
        assert "small" in method


class TestFastaParsing:
    """Test FASTA file parsing."""

    def test_parse_simple_fasta(self) -> None:
        """Test parsing simple FASTA file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">chr1\n")
            f.write("ATGCATGCATGC\n")
            f.write(">chr2\n")
            f.write("GCTAGCTAGCTA\n")
            f.flush()

            scaffolds = parse_fasta(Path(f.name))

            assert len(scaffolds) == 2
            assert scaffolds[0][0] == "chr1"
            assert scaffolds[0][1] == 12
            assert scaffolds[1][0] == "chr2"
            assert scaffolds[1][1] == 12

    def test_parse_multiline_fasta(self) -> None:
        """Test parsing FASTA with multiline sequences."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">scaffold1 description text\n")
            f.write("ATGC\n")
            f.write("ATGC\n")
            f.write("ATGC\n")
            f.flush()

            scaffolds = parse_fasta(Path(f.name))

            assert len(scaffolds) == 1
            assert scaffolds[0][0] == "scaffold1"  # Only first word
            assert scaffolds[0][1] == 12

    def test_parse_empty_file(self) -> None:
        """Test parsing empty FASTA file raises error."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.flush()

            with pytest.raises(ValueError, match="(No scaffolds|Input is empty)"):
                parse_fasta(Path(f.name))

    def test_parse_missing_file(self) -> None:
        """Test parsing missing file raises error."""
        with pytest.raises(FileNotFoundError):
            parse_fasta(Path("/nonexistent/file.fasta"))


class TestClassifyScaffolds:
    """Test full scaffold classification."""

    @pytest.fixture
    def sample_scaffolds(self) -> list[tuple[str, int, str]]:
        """Create sample scaffolds for testing."""
        return [
            ("chr1", 100_000_000, "ATGC" * 1000),
            ("chr2", 80_000_000, "GCTA" * 1000),
            ("chr3", 60_000_000, "ATGC" * 1000),
            ("Super_scaffold_4", 40_000_000, "GCTA" * 1000),
            ("scaffold_arrow_ctg1", 1_000_000, "ATGC" * 1000),
            ("scaffold_arrow_ctg2", 500_000, "GCTA" * 1000),
        ]

    def test_classify_chromosomes(
        self, sample_scaffolds: list[tuple[str, int, str]]
    ) -> None:
        """Test classification identifies chromosomes."""
        results, stats = classify_scaffolds(sample_scaffolds)

        # Should identify chr1, chr2, chr3, Super_scaffold_4 as chromosomes
        chromosomes = [r for r in results if r.classification == "chromosome"]
        assert len(chromosomes) >= 4

    def test_classify_fragments(
        self, sample_scaffolds: list[tuple[str, int, str]]
    ) -> None:
        """Test classification identifies fragments."""
        results, stats = classify_scaffolds(sample_scaffolds)

        # scaffold_arrow_ctg should be unplaced
        ctg_results = [r for r in results if "ctg" in r.name]
        assert all(r.classification == "unplaced" for r in ctg_results)

    def test_stats_calculation(
        self, sample_scaffolds: list[tuple[str, int, str]]
    ) -> None:
        """Test statistics are calculated correctly."""
        results, stats = classify_scaffolds(sample_scaffolds)

        assert stats.total_scaffolds == 6
        assert stats.total_length == sum(s[1] for s in sample_scaffolds)
        assert stats.n50 > 0
        assert stats.chromosome_count > 0

    def test_karyotype_adjustment(
        self, sample_scaffolds: list[tuple[str, int, str]]
    ) -> None:
        """Test karyotype adjustment."""
        # Force 3 chromosomes
        results, stats = classify_scaffolds(
            sample_scaffolds, expected_chromosomes=3
        )

        chromosomes = [r for r in results if r.classification == "chromosome"]
        assert len(chromosomes) == 3

    def test_empty_scaffolds_raises(self) -> None:
        """Test empty scaffold list raises error."""
        with pytest.raises(ValueError, match="No scaffolds"):
            classify_scaffolds([])


class TestScaffoldInfo:
    """Test ScaffoldInfo dataclass."""

    def test_to_dict(self) -> None:
        """Test conversion to dictionary."""
        info = ScaffoldInfo(
            name="chr1",
            length=100_000_000,
            classification="chromosome",
            confidence=0.95,
            detection_method="name_chr_explicit",
            chromosome_id="1",
        )
        d = info.to_dict()

        assert d["name"] == "chr1"
        assert d["length"] == 100_000_000
        assert d["classification"] == "chromosome"
        assert d["chromosome_id"] == "1"


class TestAssemblyStats:
    """Test AssemblyStats dataclass."""

    def test_to_dict(self) -> None:
        """Test conversion to dictionary."""
        stats = AssemblyStats(
            total_scaffolds=100,
            total_length=3_000_000_000,
            n50=50_000_000,
            n90=10_000_000,
            chromosome_count=24,
            chromosome_length=2_900_000_000,
            chromosome_n50=120_000_000,
            unlocalized_count=10,
            unplaced_count=66,
            largest_scaffold=250_000_000,
            gc_content=0.42,
        )
        d = stats.to_dict()

        assert d["total_scaffolds"] == 100
        assert d["gc_content"] == 0.42
