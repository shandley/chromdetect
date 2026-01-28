"""
Tests for multi-assembly QC dashboard functionality.
"""

import json
import tempfile
from pathlib import Path

import pytest

from chromdetect.core import AssemblyStats, ScaffoldInfo
from chromdetect.dashboard import (
    AssemblyQC,
    DashboardResult,
    _check_qc_issues,
    analyze_assembly,
    analyze_multiple_assemblies,
    format_dashboard_summary,
    format_dashboard_tsv,
    generate_dashboard_html,
)


# --- Fixtures ---


@pytest.fixture
def sample_stats() -> AssemblyStats:
    """Sample assembly statistics."""
    return AssemblyStats(
        total_scaffolds=100,
        total_length=3_000_000_000,
        n50=150_000_000,
        n90=50_000_000,
        chromosome_count=22,
        chromosome_length=2_800_000_000,
        chromosome_n50=140_000_000,
        unlocalized_count=5,
        unplaced_count=73,
        largest_scaffold=250_000_000,
    )


@pytest.fixture
def sample_scaffolds() -> list[ScaffoldInfo]:
    """Sample scaffold list."""
    scaffolds = []
    for i in range(1, 23):
        scaffolds.append(
            ScaffoldInfo(
                name=f"chr{i}",
                length=100_000_000 + i * 5_000_000,
                classification="chromosome",
                confidence=0.99,
                detection_method="name_ucsc",
                chromosome_id=str(i),
            )
        )
    return scaffolds


@pytest.fixture
def sample_qc(sample_stats, sample_scaffolds) -> AssemblyQC:
    """Sample AssemblyQC object."""
    return AssemblyQC(
        name="test_assembly",
        path="/path/to/test_assembly.fasta",
        stats=sample_stats,
        scaffolds=sample_scaffolds,
        warnings=["Low chromosome count (22)"],
        flags=[],
    )


# --- Test AssemblyQC ---


class TestAssemblyQC:
    """Tests for AssemblyQC dataclass."""

    def test_chromosome_pct(self, sample_stats, sample_scaffolds):
        """Test chromosome percentage calculation."""
        qc = AssemblyQC(
            name="test",
            path="/test",
            stats=sample_stats,
            scaffolds=sample_scaffolds,
        )

        expected = (2_800_000_000 / 3_000_000_000) * 100
        assert abs(qc.chromosome_pct - expected) < 0.01

    def test_chromosome_pct_zero_length(self, sample_scaffolds):
        """Test chromosome percentage with zero total length."""
        stats = AssemblyStats(
            total_scaffolds=0,
            total_length=0,
            n50=0,
            n90=0,
            chromosome_count=0,
            chromosome_length=0,
            chromosome_n50=0,
            unlocalized_count=0,
            unplaced_count=0,
            largest_scaffold=0,
        )
        qc = AssemblyQC(
            name="test",
            path="/test",
            stats=stats,
            scaffolds=[],
        )

        assert qc.chromosome_pct == 0.0

    def test_fragmentation_index(self, sample_stats, sample_scaffolds):
        """Test fragmentation index calculation."""
        qc = AssemblyQC(
            name="test",
            path="/test",
            stats=sample_stats,
            scaffolds=sample_scaffolds,
        )

        expected = 100 / 22  # total_scaffolds / chromosome_count
        assert abs(qc.fragmentation_index - expected) < 0.01

    def test_fragmentation_index_no_chromosomes(self):
        """Test fragmentation index with no chromosomes."""
        stats = AssemblyStats(
            total_scaffolds=100,
            total_length=1_000_000,
            n50=10_000,
            n90=1_000,
            chromosome_count=0,
            chromosome_length=0,
            chromosome_n50=0,
            unlocalized_count=0,
            unplaced_count=100,
            largest_scaffold=50_000,
        )
        qc = AssemblyQC(
            name="test",
            path="/test",
            stats=stats,
            scaffolds=[],
        )

        assert qc.fragmentation_index == float("inf")

    def test_to_dict(self, sample_qc):
        """Test conversion to dictionary."""
        d = sample_qc.to_dict()

        assert d["name"] == "test_assembly"
        assert d["path"] == "/path/to/test_assembly.fasta"
        assert "stats" in d
        assert "chromosome_pct" in d
        assert "fragmentation_index" in d
        assert "warnings" in d
        assert "flags" in d


# --- Test _check_qc_issues ---


class TestCheckQCIssues:
    """Tests for QC issue detection."""

    def test_no_chromosomes_flag(self):
        """Test flag for no chromosomes."""
        stats = AssemblyStats(
            total_scaffolds=100,
            total_length=1_000_000,
            n50=10_000,
            n90=1_000,
            chromosome_count=0,
            chromosome_length=0,
            chromosome_n50=0,
            unlocalized_count=0,
            unplaced_count=100,
            largest_scaffold=50_000,
        )

        warnings, flags = _check_qc_issues(stats, [])

        assert any("No chromosomes" in f for f in flags)

    def test_low_chromosome_count_warning(self):
        """Test warning for low chromosome count."""
        stats = AssemblyStats(
            total_scaffolds=10,
            total_length=100_000_000,
            n50=20_000_000,
            n90=10_000_000,
            chromosome_count=3,
            chromosome_length=90_000_000,
            chromosome_n50=30_000_000,
            unlocalized_count=0,
            unplaced_count=7,
            largest_scaffold=40_000_000,
        )

        warnings, flags = _check_qc_issues(stats, [])

        assert any("Low chromosome count" in w for w in warnings)

    def test_high_fragmentation_flag(self):
        """Test flag for high fragmentation."""
        stats = AssemblyStats(
            total_scaffolds=5000,
            total_length=1_000_000_000,
            n50=500_000,
            n90=100_000,
            chromosome_count=10,
            chromosome_length=500_000_000,
            chromosome_n50=50_000_000,
            unlocalized_count=0,
            unplaced_count=4990,
            largest_scaffold=100_000_000,
        )

        warnings, flags = _check_qc_issues(stats, [])

        assert any("fragmented" in f.lower() for f in flags)

    def test_low_n50_warning(self):
        """Test warning for low N50."""
        stats = AssemblyStats(
            total_scaffolds=1000,
            total_length=100_000_000,
            n50=500_000,  # 500kb - low
            n90=50_000,
            chromosome_count=10,
            chromosome_length=80_000_000,
            chromosome_n50=8_000_000,
            unlocalized_count=0,
            unplaced_count=990,
            largest_scaffold=20_000_000,
        )

        warnings, flags = _check_qc_issues(stats, [])

        assert any("Low N50" in w for w in warnings)

    def test_low_coverage_flag(self):
        """Test flag for low chromosome coverage."""
        stats = AssemblyStats(
            total_scaffolds=100,
            total_length=1_000_000_000,
            n50=10_000_000,
            n90=1_000_000,
            chromosome_count=10,
            chromosome_length=300_000_000,  # Only 30%
            chromosome_n50=30_000_000,
            unlocalized_count=0,
            unplaced_count=90,
            largest_scaffold=50_000_000,
        )

        warnings, flags = _check_qc_issues(stats, [])

        assert any("coverage" in f.lower() for f in flags)

    def test_good_assembly_no_issues(self, sample_stats, sample_scaffolds):
        """Test that good assembly has no major issues."""
        warnings, flags = _check_qc_issues(sample_stats, sample_scaffolds)

        # Should have no flags (but may have minor warnings)
        assert len(flags) == 0


# --- Test DashboardResult ---


class TestDashboardResult:
    """Tests for DashboardResult dataclass."""

    def test_summary_computation(self, sample_qc):
        """Test that summary is computed on creation."""
        result = DashboardResult(assemblies=[sample_qc])

        assert result.summary is not None
        assert result.summary["assembly_count"] == 1
        assert "n50" in result.summary
        assert "chromosome_count" in result.summary

    def test_generated_at_timestamp(self, sample_qc):
        """Test that timestamp is generated."""
        result = DashboardResult(assemblies=[sample_qc])

        assert result.generated_at != ""
        assert "T" in result.generated_at  # ISO format

    def test_get_rankings(self, sample_stats, sample_scaffolds):
        """Test ranking generation."""
        qc1 = AssemblyQC(
            name="best",
            path="/best",
            stats=sample_stats,
            scaffolds=sample_scaffolds,
        )

        low_stats = AssemblyStats(
            total_scaffolds=1000,
            total_length=1_000_000_000,
            n50=1_000_000,
            n90=100_000,
            chromosome_count=5,
            chromosome_length=500_000_000,
            chromosome_n50=100_000_000,
            unlocalized_count=0,
            unplaced_count=995,
            largest_scaffold=200_000_000,
        )
        qc2 = AssemblyQC(
            name="worst",
            path="/worst",
            stats=low_stats,
            scaffolds=[],
        )

        result = DashboardResult(assemblies=[qc1, qc2])
        rankings = result.get_rankings()

        # Best N50 should be first
        assert rankings["by_n50"][0] == "best"
        assert rankings["by_chromosome_count"][0] == "best"

    def test_to_dict(self, sample_qc):
        """Test conversion to dictionary."""
        result = DashboardResult(assemblies=[sample_qc])
        d = result.to_dict()

        assert "generated_at" in d
        assert "summary" in d
        assert "assemblies" in d
        assert len(d["assemblies"]) == 1

        # Verify JSON serializable
        json_str = json.dumps(d)
        parsed = json.loads(json_str)
        assert parsed["summary"]["assembly_count"] == 1


# --- Test analyze_assembly ---


class TestAnalyzeAssembly:
    """Tests for single assembly analysis."""

    def test_analyze_simple_assembly(self):
        """Test analyzing a simple FASTA file."""
        content = ">chr1\n" + "A" * 50_000_000 + "\n>chr2\n" + "T" * 40_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "test.fasta"
            fasta_path.write_text(content)

            qc = analyze_assembly(fasta_path, min_chromosome_size=10_000_000)

            assert qc.name == "test"
            assert qc.stats.total_scaffolds == 2
            assert qc.stats.chromosome_count == 2
            assert qc.stats.total_length == 90_000_000

    def test_analyze_with_expected_chromosomes(self):
        """Test analyzing with expected chromosome count."""
        content = ">scaffold_1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "test.fasta"
            fasta_path.write_text(content)

            qc = analyze_assembly(
                fasta_path,
                min_chromosome_size=10_000_000,
                expected_chromosomes=1,
            )

            assert qc.stats.chromosome_count >= 0


# --- Test analyze_multiple_assemblies ---


class TestAnalyzeMultipleAssemblies:
    """Tests for multi-assembly analysis."""

    def test_analyze_multiple(self):
        """Test analyzing multiple FASTA files."""
        content1 = ">chr1\n" + "A" * 50_000_000 + "\n"
        content2 = ">chr1\n" + "T" * 60_000_000 + "\n>chr2\n" + "G" * 40_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta1 = Path(tmpdir) / "assembly1.fasta"
            fasta2 = Path(tmpdir) / "assembly2.fasta"

            fasta1.write_text(content1)
            fasta2.write_text(content2)

            result = analyze_multiple_assemblies([fasta1, fasta2])

            assert len(result.assemblies) == 2
            assert result.assemblies[0].name == "assembly1"
            assert result.assemblies[1].name == "assembly2"

    def test_handles_failed_assembly(self):
        """Test that failed assemblies are handled gracefully."""
        content = ">chr1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            good_fasta = Path(tmpdir) / "good.fasta"
            bad_fasta = Path(tmpdir) / "bad.fasta"

            good_fasta.write_text(content)
            bad_fasta.write_text("not a valid fasta")

            result = analyze_multiple_assemblies([good_fasta, bad_fasta])

            assert len(result.assemblies) == 2
            # Bad assembly should have a flag
            bad_qc = next(a for a in result.assemblies if a.name == "bad")
            assert any("Failed" in f for f in bad_qc.flags)


# --- Test formatters ---


class TestFormatters:
    """Tests for output formatters."""

    def test_format_dashboard_summary(self, sample_qc):
        """Test human-readable summary formatting."""
        result = DashboardResult(assemblies=[sample_qc])
        output = format_dashboard_summary(result)

        assert "MULTI-ASSEMBLY QC DASHBOARD" in output
        assert "test_assembly" in output
        assert "N50" in output
        assert "Chromosomes" in output

    def test_format_dashboard_tsv(self, sample_qc):
        """Test TSV formatting."""
        result = DashboardResult(assemblies=[sample_qc])
        output = format_dashboard_tsv(result)
        lines = output.strip().split("\n")

        # Check header
        header = lines[0].split("\t")
        assert "name" in header
        assert "n50" in header
        assert "chromosomes" in header

        # Check data row
        assert len(lines) == 2
        data = lines[1].split("\t")
        assert data[0] == "test_assembly"

    def test_generate_dashboard_html(self, sample_qc):
        """Test HTML dashboard generation."""
        result = DashboardResult(assemblies=[sample_qc])
        output = generate_dashboard_html(result)

        assert "<!DOCTYPE html>" in output
        assert "Multi-Assembly QC Dashboard" in output
        assert "test_assembly" in output
        assert "chart.js" in output.lower()
        assert "n50Chart" in output


# --- Test CLI integration ---


class TestCLIIntegration:
    """Integration tests for CLI dashboard mode."""

    def test_cli_dashboard_summary(self, capsys, monkeypatch):
        """Test CLI with --dashboard flag."""
        content1 = ">chr1\n" + "A" * 50_000_000 + "\n"
        content2 = ">chr1\n" + "T" * 60_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta1 = Path(tmpdir) / "a1.fasta"
            fasta2 = Path(tmpdir) / "a2.fasta"

            fasta1.write_text(content1)
            fasta2.write_text(content2)

            import sys
            from chromdetect.cli import main

            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", "--dashboard", str(fasta1), str(fasta2), "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            assert exc_info.value.code == 0

            captured = capsys.readouterr()
            assert "MULTI-ASSEMBLY QC DASHBOARD" in captured.out

    def test_cli_dashboard_json(self, capsys, monkeypatch):
        """Test CLI dashboard with JSON output."""
        content = ">chr1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "test.fasta"
            fasta.write_text(content)

            import sys
            from chromdetect.cli import main

            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", "--dashboard", str(fasta), "--format", "json", "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            assert exc_info.value.code == 0

            captured = capsys.readouterr()
            parsed = json.loads(captured.out)
            assert "assemblies" in parsed
            assert "summary" in parsed

    def test_cli_dashboard_html(self, capsys, monkeypatch):
        """Test CLI dashboard with HTML output."""
        content = ">chr1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "test.fasta"
            fasta.write_text(content)

            import sys
            from chromdetect.cli import main

            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", "--dashboard", str(fasta), "--format", "html", "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            assert exc_info.value.code == 0

            captured = capsys.readouterr()
            assert "<!DOCTYPE html>" in captured.out
            assert "chart.js" in captured.out.lower()

    def test_cli_dashboard_to_file(self, monkeypatch):
        """Test CLI dashboard output to file."""
        content = ">chr1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "test.fasta"
            output = Path(tmpdir) / "dashboard.html"

            fasta.write_text(content)

            import sys
            from chromdetect.cli import main

            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", "--dashboard", str(fasta), "-o", str(output),
                 "--format", "html", "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            assert exc_info.value.code == 0
            assert output.exists()

            html_content = output.read_text()
            assert "<!DOCTYPE html>" in html_content

    def test_cli_dashboard_missing_file(self, capsys, monkeypatch):
        """Test CLI error handling for missing file."""
        import sys
        from chromdetect.cli import main

        monkeypatch.setattr(
            sys, "argv",
            ["chromdetect", "--dashboard", "/nonexistent/file.fasta", "-q"]
        )

        with pytest.raises(SystemExit) as exc_info:
            main()

        assert exc_info.value.code == 66  # EXIT_NOINPUT

    def test_cli_dashboard_tsv(self, capsys, monkeypatch):
        """Test CLI dashboard with TSV output."""
        content = ">chr1\n" + "A" * 50_000_000 + "\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "test.fasta"
            fasta.write_text(content)

            import sys
            from chromdetect.cli import main

            monkeypatch.setattr(
                sys, "argv",
                ["chromdetect", "--dashboard", str(fasta), "--format", "tsv", "-q"]
            )

            with pytest.raises(SystemExit) as exc_info:
                main()

            assert exc_info.value.code == 0

            captured = capsys.readouterr()
            lines = captured.out.strip().split("\n")
            assert "name" in lines[0]
            assert "n50" in lines[0]
