"""
Tests for HTML report generation functionality.
"""



from chromdetect.core import AssemblyStats, ScaffoldInfo
from chromdetect.html_report import (
    _format_bp,
    _generate_bar_chart,
    _generate_pie_chart,
    generate_html_report,
)


def make_scaffold(name, length, classification, confidence=0.9, chr_id=None):
    """Helper to create ScaffoldInfo for testing."""
    return ScaffoldInfo(
        name=name,
        length=length,
        classification=classification,
        confidence=confidence,
        detection_method="test",
        chromosome_id=chr_id,
        gc_content=0.4,
    )


def make_stats():
    """Helper to create AssemblyStats for testing."""
    return AssemblyStats(
        total_scaffolds=100,
        total_length=2_800_000_000,
        n50=150_000_000,
        n90=80_000_000,
        chromosome_count=23,
        chromosome_length=2_700_000_000,
        chromosome_n50=180_000_000,
        unlocalized_count=5,
        unplaced_count=72,
        largest_scaffold=250_000_000,
        gc_content=0.41,
    )


class TestFormatBP:
    """Tests for _format_bp function."""

    def test_gigabases(self):
        """Test formatting gigabases."""
        assert _format_bp(2_800_000_000) == "2.80 Gb"
        assert _format_bp(1_000_000_000) == "1.00 Gb"

    def test_megabases(self):
        """Test formatting megabases."""
        assert _format_bp(150_000_000) == "150.00 Mb"
        assert _format_bp(1_000_000) == "1.00 Mb"

    def test_kilobases(self):
        """Test formatting kilobases."""
        assert _format_bp(50_000) == "50.00 kb"
        assert _format_bp(1_000) == "1.00 kb"

    def test_bases(self):
        """Test formatting bases."""
        assert _format_bp(500) == "500 bp"
        assert _format_bp(1) == "1 bp"


class TestGeneratePieChart:
    """Tests for pie chart generation."""

    def test_empty_data(self):
        """Test pie chart with no data."""
        result = _generate_pie_chart([], "Test Chart")
        assert "No data" in result

    def test_all_zeros(self):
        """Test pie chart with all zero values."""
        data = [("A", 0, "#ff0000"), ("B", 0, "#00ff00")]
        result = _generate_pie_chart(data, "Test Chart")
        assert "No data" in result

    def test_valid_data(self):
        """Test pie chart with valid data."""
        data = [
            ("Chromosomes", 23, "#4CAF50"),
            ("Unplaced", 77, "#9E9E9E"),
        ]
        result = _generate_pie_chart(data, "Scaffold Types")

        assert "<svg" in result
        assert "</svg>" in result
        assert "Scaffold Types" in result

    def test_single_slice(self):
        """Test pie chart with single value (full circle)."""
        data = [("Only", 100, "#ff0000")]
        result = _generate_pie_chart(data, "Single")

        assert "<svg" in result
        # Should draw a circle for 100%
        assert "<circle" in result or "<path" in result


class TestGenerateBarChart:
    """Tests for bar chart generation."""

    def test_empty_data(self):
        """Test bar chart with no data."""
        result = _generate_bar_chart([], "Test Chart")
        assert "No data" in result

    def test_valid_data(self):
        """Test bar chart with valid data."""
        data = [
            ("chr1", 200_000_000, "#4CAF50"),
            ("chr2", 180_000_000, "#4CAF50"),
            ("chr3", 160_000_000, "#4CAF50"),
        ]
        result = _generate_bar_chart(data, "Chromosome Sizes")

        assert "<svg" in result
        assert "</svg>" in result
        assert "Chromosome Sizes" in result
        assert "<rect" in result  # Bars

    def test_max_bars_limit(self):
        """Test that max_bars limits the number of bars."""
        data = [(f"chr{i}", i * 1000, "#4CAF50") for i in range(50)]
        result = _generate_bar_chart(data, "Many Bars", max_bars=10)

        # Should only include first 10 bars
        assert "chr0" in result
        assert "chr9" in result
        # chr10 onwards should not be included
        assert "chr10" not in result


class TestGenerateHTMLReport:
    """Tests for full HTML report generation."""

    def test_basic_report(self):
        """Test basic HTML report generation."""
        results = [
            make_scaffold("chr1", 200_000_000, "chromosome", chr_id="1"),
            make_scaffold("chr2", 180_000_000, "chromosome", chr_id="2"),
            make_scaffold("scaffold1", 50_000, "unplaced"),
        ]
        stats = make_stats()

        html = generate_html_report(results, stats, "Test Assembly")

        assert "<!DOCTYPE html>" in html
        assert "<html" in html
        assert "</html>" in html
        assert "Test Assembly" in html
        assert "ChromDetect" in html

    def test_report_contains_stats(self):
        """Test that report contains statistics."""
        results = [make_scaffold("chr1", 100_000_000, "chromosome")]
        stats = make_stats()

        html = generate_html_report(results, stats, "Test")

        assert "2.80 Gb" in html or "2,800,000,000" in html
        assert "23" in html  # chromosome count
        assert "0.85" in html or "85" in html  # quality score

    def test_report_contains_charts(self):
        """Test that report contains SVG charts."""
        results = [
            make_scaffold("chr1", 200_000_000, "chromosome"),
            make_scaffold("scaffold1", 50_000, "unplaced"),
        ]
        stats = make_stats()

        html = generate_html_report(results, stats, "Test")

        # Should have SVG charts
        assert "<svg" in html
        assert "</svg>" in html

    def test_report_contains_table(self):
        """Test that report contains scaffold table."""
        results = [
            make_scaffold("chr1", 200_000_000, "chromosome", chr_id="1"),
            make_scaffold("chr2", 180_000_000, "chromosome", chr_id="2"),
        ]
        stats = make_stats()

        html = generate_html_report(results, stats, "Test")

        assert "<table" in html
        assert "</table>" in html
        assert "chr1" in html
        assert "chr2" in html

    def test_report_has_styling(self):
        """Test that report has CSS styling."""
        results = [make_scaffold("chr1", 100_000_000, "chromosome")]
        stats = make_stats()

        html = generate_html_report(results, stats, "Test")

        assert "<style>" in html
        assert "</style>" in html

    def test_report_html_escaping(self):
        """Test that special characters are escaped."""
        results = [
            make_scaffold("scaffold<test>&special", 100_000, "unplaced"),
        ]
        stats = make_stats()

        html = generate_html_report(results, stats, "Test<>")

        # Should be escaped
        assert "&lt;" in html or "scaffold&lt;test&gt;" in html
        assert "<script>" not in html  # No XSS

    def test_report_limits_scaffolds(self):
        """Test that table limits to top 100 scaffolds."""
        results = [
            make_scaffold(f"scaffold{i}", 1000 + i, "unplaced")
            for i in range(200)
        ]
        stats = make_stats()

        html = generate_html_report(results, stats, "Test")

        # Should only include up to 100 scaffolds
        # (sorted by size descending)
        assert "scaffold199" in html  # Largest should be included
        assert html.count("<tr>") <= 101  # Header + 100 rows
