"""
Tests for assembly version tracking functionality.
"""

import json
import tempfile
from pathlib import Path

import pytest

from chromdetect.core import AssemblyStats, ScaffoldInfo
from chromdetect.version import (
    ChangeType,
    ScaffoldChange,
    VersionComparisonResult,
    _detect_change_type,
    _match_scaffolds,
    compare_assemblies,
    compare_fasta_files,
    format_version_json,
    format_version_summary,
    format_version_tsv,
)


# --- Fixtures ---


@pytest.fixture
def sample_scaffolds_v1() -> list[ScaffoldInfo]:
    """Sample scaffolds for version 1."""
    return [
        ScaffoldInfo(
            name="chr1",
            length=100_000_000,
            classification="chromosome",
            confidence=0.99,
            detection_method="name_ucsc",
            chromosome_id="1",
        ),
        ScaffoldInfo(
            name="chr2",
            length=80_000_000,
            classification="chromosome",
            confidence=0.99,
            detection_method="name_ucsc",
            chromosome_id="2",
        ),
        ScaffoldInfo(
            name="scaffold_100",
            length=5_000_000,
            classification="unplaced",
            confidence=0.6,
            detection_method="size_small",
        ),
        ScaffoldInfo(
            name="scaffold_101",
            length=3_000_000,
            classification="unplaced",
            confidence=0.6,
            detection_method="size_small",
        ),
    ]


@pytest.fixture
def sample_stats_v1() -> AssemblyStats:
    """Sample stats for version 1."""
    return AssemblyStats(
        total_scaffolds=4,
        total_length=188_000_000,
        n50=80_000_000,
        n90=5_000_000,
        chromosome_count=2,
        chromosome_length=180_000_000,
        chromosome_n50=100_000_000,
        unlocalized_count=0,
        unplaced_count=2,
        largest_scaffold=100_000_000,
    )


@pytest.fixture
def sample_scaffolds_v2() -> list[ScaffoldInfo]:
    """Sample scaffolds for version 2 with some changes."""
    return [
        ScaffoldInfo(
            name="chr1",
            length=102_000_000,  # Slightly larger (2% growth)
            classification="chromosome",
            confidence=0.99,
            detection_method="name_ucsc",
            chromosome_id="1",
        ),
        ScaffoldInfo(
            name="chr2",
            length=80_000_000,  # Same
            classification="chromosome",
            confidence=0.99,
            detection_method="name_ucsc",
            chromosome_id="2",
        ),
        ScaffoldInfo(
            name="chr3",  # Promoted from scaffold_100
            length=5_000_000,
            classification="chromosome",
            confidence=0.8,
            detection_method="name_ucsc",
            chromosome_id="3",
        ),
        # scaffold_101 is removed
        ScaffoldInfo(
            name="scaffold_200",  # New scaffold
            length=2_000_000,
            classification="unplaced",
            confidence=0.6,
            detection_method="size_small",
        ),
    ]


@pytest.fixture
def sample_stats_v2() -> AssemblyStats:
    """Sample stats for version 2."""
    return AssemblyStats(
        total_scaffolds=4,
        total_length=189_000_000,
        n50=80_000_000,
        n90=5_000_000,
        chromosome_count=3,
        chromosome_length=187_000_000,
        chromosome_n50=80_000_000,
        unlocalized_count=0,
        unplaced_count=1,
        largest_scaffold=102_000_000,
    )


# --- Test _match_scaffolds ---


class TestMatchScaffolds:
    """Tests for scaffold matching between versions."""

    def test_exact_name_match(self):
        """Scaffolds with same names should match exactly."""
        v1 = [
            ScaffoldInfo("chr1", 100, "chromosome", 0.9, "test"),
            ScaffoldInfo("chr2", 200, "chromosome", 0.9, "test"),
        ]
        v2 = [
            ScaffoldInfo("chr1", 100, "chromosome", 0.9, "test"),
            ScaffoldInfo("chr2", 200, "chromosome", 0.9, "test"),
        ]

        mapping, removed, added = _match_scaffolds(v1, v2)

        assert mapping == {"chr1": "chr1", "chr2": "chr2"}
        assert removed == set()
        assert added == set()

    def test_detect_removed_scaffolds(self):
        """Scaffolds only in v1 should be detected as removed."""
        v1 = [
            ScaffoldInfo("chr1", 100, "chromosome", 0.9, "test"),
            ScaffoldInfo("removed_scaffold", 50, "unplaced", 0.6, "test"),
        ]
        v2 = [
            ScaffoldInfo("chr1", 100, "chromosome", 0.9, "test"),
        ]

        mapping, removed, added = _match_scaffolds(v1, v2)

        assert "chr1" in mapping
        assert "removed_scaffold" in removed
        assert added == set()

    def test_detect_added_scaffolds(self):
        """Scaffolds only in v2 should be detected as added."""
        v1 = [
            ScaffoldInfo("chr1", 100, "chromosome", 0.9, "test"),
        ]
        v2 = [
            ScaffoldInfo("chr1", 100, "chromosome", 0.9, "test"),
            ScaffoldInfo("new_scaffold", 50, "unplaced", 0.6, "test"),
        ]

        mapping, removed, added = _match_scaffolds(v1, v2)

        assert "chr1" in mapping
        assert removed == set()
        assert "new_scaffold" in added

    def test_match_by_chromosome_id(self):
        """Scaffolds with same chromosome ID should match even if renamed."""
        v1 = [
            ScaffoldInfo("scaffold_1", 100_000_000, "chromosome", 0.9, "test", chromosome_id="1"),
        ]
        v2 = [
            ScaffoldInfo("chr1", 100_000_000, "chromosome", 0.9, "test", chromosome_id="1"),
        ]

        mapping, removed, added = _match_scaffolds(v1, v2)

        assert mapping == {"scaffold_1": "chr1"}
        assert removed == set()
        assert added == set()

    def test_match_by_size_within_tolerance(self):
        """Scaffolds with similar size and classification should match."""
        v1 = [
            ScaffoldInfo("old_name", 100_000_000, "chromosome", 0.9, "test", chromosome_id="1"),
        ]
        v2 = [
            ScaffoldInfo("new_name", 105_000_000, "chromosome", 0.9, "test", chromosome_id="1"),  # 5% larger
        ]

        # Default tolerance is 10%
        mapping, removed, added = _match_scaffolds(v1, v2, size_tolerance=0.1)

        assert mapping == {"old_name": "new_name"}

    def test_no_match_beyond_tolerance(self):
        """Scaffolds with very different sizes should not match."""
        v1 = [
            ScaffoldInfo("scaffold_1", 100_000_000, "chromosome", 0.9, "test"),
        ]
        v2 = [
            ScaffoldInfo("scaffold_2", 50_000_000, "chromosome", 0.9, "test"),  # 50% smaller
        ]

        mapping, removed, added = _match_scaffolds(v1, v2, size_tolerance=0.1)

        # Different names, too different sizes, no chromosome_id - should not match
        assert "scaffold_1" in removed
        assert "scaffold_2" in added


# --- Test _detect_change_type ---


class TestDetectChangeType:
    """Tests for change type detection."""

    def test_added_scaffold(self):
        """New scaffold should be detected as added."""
        s2 = ScaffoldInfo("new_scaffold", 100, "unplaced", 0.6, "test")
        change_type, notes = _detect_change_type(None, s2)

        assert change_type == ChangeType.ADDED
        assert "New scaffold" in notes

    def test_removed_scaffold(self):
        """Missing scaffold should be detected as removed."""
        s1 = ScaffoldInfo("old_scaffold", 100, "unplaced", 0.6, "test")
        change_type, notes = _detect_change_type(s1, None)

        assert change_type == ChangeType.REMOVED
        assert "Removed" in notes

    def test_promoted_scaffold(self):
        """Scaffold promoted to higher classification should be detected."""
        s1 = ScaffoldInfo("scaffold_1", 100, "unplaced", 0.6, "test")
        s2 = ScaffoldInfo("scaffold_1", 100, "chromosome", 0.9, "test")

        change_type, notes = _detect_change_type(s1, s2)

        assert change_type == ChangeType.PROMOTED
        assert "unplaced" in notes
        assert "chromosome" in notes

    def test_demoted_scaffold(self):
        """Scaffold demoted to lower classification should be detected."""
        s1 = ScaffoldInfo("chr1", 100, "chromosome", 0.9, "test")
        s2 = ScaffoldInfo("chr1", 100, "unplaced", 0.6, "test")

        change_type, notes = _detect_change_type(s1, s2)

        assert change_type == ChangeType.DEMOTED
        assert "chromosome" in notes
        assert "unplaced" in notes

    def test_resized_scaffold(self):
        """Scaffold with significant size change should be detected."""
        s1 = ScaffoldInfo("chr1", 100_000_000, "chromosome", 0.9, "test")
        s2 = ScaffoldInfo("chr1", 130_000_000, "chromosome", 0.9, "test")  # 30% larger

        change_type, notes = _detect_change_type(s1, s2, size_threshold=0.2)

        assert change_type == ChangeType.RESIZED
        assert "grew" in notes

    def test_renamed_scaffold(self):
        """Scaffold with different name should be marked as renamed."""
        s1 = ScaffoldInfo("old_name", 100, "chromosome", 0.9, "test")
        s2 = ScaffoldInfo("new_name", 100, "chromosome", 0.9, "test")

        change_type, notes = _detect_change_type(s1, s2, is_renamed=True)

        assert change_type == ChangeType.RENAMED
        assert "renamed" in notes

    def test_unchanged_scaffold(self):
        """Scaffold with no changes should be marked as unchanged."""
        s1 = ScaffoldInfo("chr1", 100_000_000, "chromosome", 0.9, "test")
        s2 = ScaffoldInfo("chr1", 100_000_000, "chromosome", 0.9, "test")

        change_type, notes = _detect_change_type(s1, s2)

        assert change_type == ChangeType.UNCHANGED
        assert notes == ""

    def test_small_size_change_not_resized(self):
        """Small size changes should not be flagged as resized."""
        s1 = ScaffoldInfo("chr1", 100_000_000, "chromosome", 0.9, "test")
        s2 = ScaffoldInfo("chr1", 105_000_000, "chromosome", 0.9, "test")  # 5% larger

        change_type, notes = _detect_change_type(s1, s2, size_threshold=0.2)

        assert change_type == ChangeType.UNCHANGED


# --- Test compare_assemblies ---


class TestCompareAssemblies:
    """Tests for full assembly comparison."""

    def test_basic_comparison(
        self,
        sample_scaffolds_v1,
        sample_stats_v1,
        sample_scaffolds_v2,
        sample_stats_v2,
    ):
        """Test basic assembly comparison."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v2,
            sample_stats_v2,
            label_v1="v1.0",
            label_v2="v2.0",
        )

        assert result.version1_label == "v1.0"
        assert result.version2_label == "v2.0"
        assert isinstance(result.changes, list)
        assert len(result.changes) > 0

    def test_counts_match_changes(
        self,
        sample_scaffolds_v1,
        sample_stats_v1,
        sample_scaffolds_v2,
        sample_stats_v2,
    ):
        """Test that change counts match the actual changes."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v2,
            sample_stats_v2,
        )

        # Count changes by type
        promoted = [c for c in result.changes if c.change_type == ChangeType.PROMOTED]
        demoted = [c for c in result.changes if c.change_type == ChangeType.DEMOTED]
        added = [c for c in result.changes if c.change_type == ChangeType.ADDED]
        removed = [c for c in result.changes if c.change_type == ChangeType.REMOVED]
        unchanged = [c for c in result.changes if c.change_type == ChangeType.UNCHANGED]

        assert result.promoted_count == len(promoted)
        assert result.demoted_count == len(demoted)
        assert result.added_count == len(added)
        assert result.removed_count == len(removed)
        assert result.unchanged_count == len(unchanged)

    def test_n50_change_calculated(
        self,
        sample_scaffolds_v1,
        sample_stats_v1,
        sample_scaffolds_v2,
        sample_stats_v2,
    ):
        """Test that N50 change is correctly calculated."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v2,
            sample_stats_v2,
        )

        expected_n50_change = sample_stats_v2.n50 - sample_stats_v1.n50
        assert result.n50_change == expected_n50_change

    def test_chromosome_count_change(
        self,
        sample_scaffolds_v1,
        sample_stats_v1,
        sample_scaffolds_v2,
        sample_stats_v2,
    ):
        """Test that chromosome count change is correctly calculated."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v2,
            sample_stats_v2,
        )

        expected_chr_change = sample_stats_v2.chromosome_count - sample_stats_v1.chromosome_count
        assert result.chromosome_count_change == expected_chr_change

    def test_summary_generated(
        self,
        sample_scaffolds_v1,
        sample_stats_v1,
        sample_scaffolds_v2,
        sample_stats_v2,
    ):
        """Test that a summary is generated."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v2,
            sample_stats_v2,
        )

        assert result.summary != ""
        assert isinstance(result.summary, str)

    def test_identical_assemblies(self, sample_scaffolds_v1, sample_stats_v1):
        """Test comparison of identical assemblies."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v1,  # Same scaffolds
            sample_stats_v1,  # Same stats
        )

        # All scaffolds should be unchanged
        assert result.unchanged_count == len(sample_scaffolds_v1)
        assert result.promoted_count == 0
        assert result.demoted_count == 0
        assert result.added_count == 0
        assert result.removed_count == 0


# --- Test compare_fasta_files ---


class TestCompareFastaFiles:
    """Tests for FASTA file comparison."""

    def test_compare_fasta_files(self):
        """Test comparing two FASTA files."""
        v1_content = ">chr1\n" + "A" * 50_000_000 + "\n>chr2\n" + "T" * 30_000_000 + "\n"
        v2_content = ">chr1\n" + "A" * 52_000_000 + "\n>chr2\n" + "T" * 30_000_000 + "\n>chr3\n" + "G" * 20_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            v1_path = Path(tmpdir) / "v1.fasta"
            v2_path = Path(tmpdir) / "v2.fasta"

            v1_path.write_text(v1_content)
            v2_path.write_text(v2_content)

            result = compare_fasta_files(
                v1_path,
                v2_path,
                min_chromosome_size=10_000_000,
            )

            assert result.version1_label == "v1.fasta"
            assert result.version2_label == "v2.fasta"
            assert result.added_count >= 1  # chr3 is added
            assert len(result.changes) == 3

    def test_compare_fasta_with_custom_labels(self):
        """Test comparing FASTA files with custom labels."""
        content = ">chr1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            v1_path = Path(tmpdir) / "assembly.fasta"
            v2_path = Path(tmpdir) / "assembly_updated.fasta"

            v1_path.write_text(content)
            v2_path.write_text(content)

            result = compare_fasta_files(
                v1_path,
                v2_path,
                label_v1="Assembly v1.0",
                label_v2="Assembly v2.0",
            )

            assert result.version1_label == "Assembly v1.0"
            assert result.version2_label == "Assembly v2.0"


# --- Test output formatters ---


class TestFormatters:
    """Tests for output formatters."""

    def test_format_version_summary(
        self,
        sample_scaffolds_v1,
        sample_stats_v1,
        sample_scaffolds_v2,
        sample_stats_v2,
    ):
        """Test human-readable summary formatting."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v2,
            sample_stats_v2,
            label_v1="v1.0",
            label_v2="v2.0",
        )

        output = format_version_summary(result)

        assert "ASSEMBLY VERSION COMPARISON" in output
        assert "v1.0" in output
        assert "v2.0" in output
        assert "METRIC CHANGES" in output
        assert "SCAFFOLD CHANGES" in output
        assert "N50" in output

    def test_format_version_tsv(
        self,
        sample_scaffolds_v1,
        sample_stats_v1,
        sample_scaffolds_v2,
        sample_stats_v2,
    ):
        """Test TSV formatting."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v2,
            sample_stats_v2,
        )

        output = format_version_tsv(result)
        lines = output.strip().split("\n")

        # Check header
        header = lines[0].split("\t")
        assert "change_type" in header
        assert "name_v1" in header
        assert "name_v2" in header
        assert "class_v1" in header
        assert "class_v2" in header

        # Check data rows
        assert len(lines) > 1  # Header + at least one data row

    def test_format_version_json(
        self,
        sample_scaffolds_v1,
        sample_stats_v1,
        sample_scaffolds_v2,
        sample_stats_v2,
    ):
        """Test JSON formatting."""
        result = compare_assemblies(
            sample_scaffolds_v1,
            sample_stats_v1,
            sample_scaffolds_v2,
            sample_stats_v2,
            label_v1="v1.0",
            label_v2="v2.0",
        )

        output = format_version_json(result)

        assert isinstance(output, dict)
        assert output["version1"] == "v1.0"
        assert output["version2"] == "v2.0"
        assert "stats_v1" in output
        assert "stats_v2" in output
        assert "metrics" in output
        assert "counts" in output
        assert "changes" in output
        assert "summary" in output

        # Verify JSON is serializable
        json_str = json.dumps(output)
        parsed = json.loads(json_str)
        assert parsed["version1"] == "v1.0"


# --- Test ScaffoldChange ---


class TestScaffoldChange:
    """Tests for ScaffoldChange dataclass."""

    def test_primary_name_prefers_v2(self):
        """Primary name should prefer v2 name."""
        change = ScaffoldChange(
            name_v1="old_name",
            name_v2="new_name",
            change_type=ChangeType.RENAMED,
            class_v1="chromosome",
            class_v2="chromosome",
            length_v1=100,
            length_v2=100,
        )

        assert change.primary_name == "new_name"

    def test_primary_name_falls_back_to_v1(self):
        """Primary name should fall back to v1 if v2 is None."""
        change = ScaffoldChange(
            name_v1="removed_scaffold",
            name_v2=None,
            change_type=ChangeType.REMOVED,
            class_v1="unplaced",
            class_v2=None,
            length_v1=100,
            length_v2=0,
        )

        assert change.primary_name == "removed_scaffold"

    def test_primary_name_handles_both_none(self):
        """Primary name should return 'unknown' if both are None."""
        change = ScaffoldChange(
            name_v1=None,
            name_v2=None,
            change_type=ChangeType.UNCHANGED,
            class_v1=None,
            class_v2=None,
            length_v1=0,
            length_v2=0,
        )

        assert change.primary_name == "unknown"


# --- Test VersionComparisonResult ---


class TestVersionComparisonResult:
    """Tests for VersionComparisonResult dataclass."""

    def test_n50_change_pct(self, sample_stats_v1, sample_stats_v2):
        """Test N50 percent change calculation."""
        result = VersionComparisonResult(
            version1_label="v1",
            version2_label="v2",
            stats_v1=sample_stats_v1,
            stats_v2=sample_stats_v2,
        )

        expected_pct = (result.n50_change / sample_stats_v1.n50) * 100
        assert result.n50_change_pct == expected_pct

    def test_n50_change_pct_zero_division(self, sample_stats_v2):
        """Test N50 percent change with zero v1 N50."""
        zero_n50_stats = AssemblyStats(
            total_scaffolds=1,
            total_length=100,
            n50=0,
            n90=0,
            chromosome_count=0,
            chromosome_length=0,
            chromosome_n50=0,
            unlocalized_count=0,
            unplaced_count=1,
            largest_scaffold=100,
        )

        result = VersionComparisonResult(
            version1_label="v1",
            version2_label="v2",
            stats_v1=zero_n50_stats,
            stats_v2=sample_stats_v2,
        )

        assert result.n50_change_pct == 0.0


# --- Test CLI integration ---


class TestCLIIntegration:
    """Integration tests for CLI version comparison."""

    def test_cli_compare_versions(self, capsys, monkeypatch):
        """Test CLI with --compare-versions flag."""
        v1_content = ">chr1\n" + "A" * 50_000_000 + "\n"
        v2_content = ">chr1\n" + "A" * 52_000_000 + "\n>chr2\n" + "T" * 30_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            v1_path = Path(tmpdir) / "v1.fasta"
            v2_path = Path(tmpdir) / "v2.fasta"

            v1_path.write_text(v1_content)
            v2_path.write_text(v2_content)

            import sys
            from chromdetect.cli import main

            # Test summary format
            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", str(v1_path), "--compare-versions", str(v2_path), "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            assert exc_info.value.code == 0

            captured = capsys.readouterr()
            assert "ASSEMBLY VERSION COMPARISON" in captured.out

    def test_cli_compare_versions_json(self, capsys, monkeypatch):
        """Test CLI version comparison with JSON output."""
        v1_content = ">chr1\n" + "A" * 50_000_000 + "\n"
        v2_content = ">chr1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            v1_path = Path(tmpdir) / "v1.fasta"
            v2_path = Path(tmpdir) / "v2.fasta"

            v1_path.write_text(v1_content)
            v2_path.write_text(v2_content)

            import sys
            from chromdetect.cli import main

            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", str(v1_path), "--compare-versions", str(v2_path),
                 "--format", "json", "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            assert exc_info.value.code == 0

            captured = capsys.readouterr()
            parsed = json.loads(captured.out)
            assert "version1" in parsed
            assert "version2" in parsed

    def test_cli_compare_versions_tsv(self, capsys, monkeypatch):
        """Test CLI version comparison with TSV output."""
        v1_content = ">chr1\n" + "A" * 50_000_000 + "\n"
        v2_content = ">chr1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            v1_path = Path(tmpdir) / "v1.fasta"
            v2_path = Path(tmpdir) / "v2.fasta"

            v1_path.write_text(v1_content)
            v2_path.write_text(v2_content)

            import sys
            from chromdetect.cli import main

            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", str(v1_path), "--compare-versions", str(v2_path),
                 "--format", "tsv", "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            assert exc_info.value.code == 0

            captured = capsys.readouterr()
            lines = captured.out.strip().split("\n")
            assert "change_type" in lines[0]

    def test_cli_compare_versions_missing_file(self, capsys, monkeypatch):
        """Test CLI error handling for missing file."""
        v1_content = ">chr1\n" + "A" * 100 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            v1_path = Path(tmpdir) / "v1.fasta"
            v2_path = Path(tmpdir) / "nonexistent.fasta"

            v1_path.write_text(v1_content)

            import sys
            from chromdetect.cli import main

            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", str(v1_path), "--compare-versions", str(v2_path), "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            # Should exit with error code 66 (NOINPUT)
            assert exc_info.value.code == 66
