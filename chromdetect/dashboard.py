"""
Multi-Assembly QC Dashboard for ChromDetect.

This module provides functionality to analyze multiple genome assemblies
and generate a comparative quality control dashboard.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from chromdetect.core import (
    AssemblyStats,
    ScaffoldInfo,
    classify_scaffolds,
    parse_fasta,
)


@dataclass
class AssemblyQC:
    """Quality control metrics for a single assembly.

    Attributes:
        name: Assembly name (typically filename)
        path: Path to the FASTA file
        stats: Assembly statistics from classification
        scaffolds: Classified scaffold information
        warnings: List of QC warnings
        flags: List of QC flags (issues that may indicate problems)
    """

    name: str
    path: str
    stats: AssemblyStats
    scaffolds: list[ScaffoldInfo]
    warnings: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)

    @property
    def chromosome_pct(self) -> float:
        """Percentage of assembly in chromosome-level scaffolds."""
        if self.stats.total_length == 0:
            return 0.0
        return (self.stats.chromosome_length / self.stats.total_length) * 100

    @property
    def fragmentation_index(self) -> float:
        """Fragmentation index (lower is better).

        Calculated as: total_scaffolds / chromosome_count
        A highly contiguous assembly has a value close to 1.
        """
        if self.stats.chromosome_count == 0:
            return float("inf")
        return self.stats.total_scaffolds / self.stats.chromosome_count

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "path": self.path,
            "stats": self.stats.to_dict(),
            "chromosome_pct": round(self.chromosome_pct, 2),
            "fragmentation_index": round(self.fragmentation_index, 2)
            if self.fragmentation_index != float("inf")
            else None,
            "warnings": self.warnings,
            "flags": self.flags,
            "scaffold_count": len(self.scaffolds),
        }


@dataclass
class DashboardResult:
    """Results from multi-assembly QC analysis.

    Attributes:
        assemblies: List of AssemblyQC objects for each assembly
        generated_at: Timestamp when dashboard was generated
        summary: Overall summary statistics
    """

    assemblies: list[AssemblyQC]
    generated_at: str = ""
    summary: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.generated_at:
            self.generated_at = datetime.now().isoformat()
        if not self.summary:
            self.summary = self._compute_summary()

    def _compute_summary(self) -> dict:
        """Compute summary statistics across all assemblies."""
        if not self.assemblies:
            return {}

        n50s = [a.stats.n50 for a in self.assemblies]
        chr_counts = [a.stats.chromosome_count for a in self.assemblies]
        total_lengths = [a.stats.total_length for a in self.assemblies]
        chr_pcts = [a.chromosome_pct for a in self.assemblies]

        return {
            "assembly_count": len(self.assemblies),
            "n50": {
                "min": min(n50s),
                "max": max(n50s),
                "mean": sum(n50s) // len(n50s),
                "best": self.assemblies[n50s.index(max(n50s))].name,
            },
            "chromosome_count": {
                "min": min(chr_counts),
                "max": max(chr_counts),
                "mean": sum(chr_counts) // len(chr_counts),
            },
            "total_length": {
                "min": min(total_lengths),
                "max": max(total_lengths),
                "mean": sum(total_lengths) // len(total_lengths),
            },
            "chromosome_pct": {
                "min": round(min(chr_pcts), 2),
                "max": round(max(chr_pcts), 2),
                "mean": round(sum(chr_pcts) / len(chr_pcts), 2),
            },
            "total_warnings": sum(len(a.warnings) for a in self.assemblies),
            "total_flags": sum(len(a.flags) for a in self.assemblies),
        }

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "generated_at": self.generated_at,
            "summary": self.summary,
            "assemblies": [a.to_dict() for a in self.assemblies],
        }

    def get_rankings(self) -> dict[str, list[str]]:
        """Get rankings of assemblies by various metrics."""
        if not self.assemblies:
            return {}

        # Sort by N50 (descending - higher is better)
        by_n50 = sorted(self.assemblies, key=lambda a: a.stats.n50, reverse=True)

        # Sort by chromosome count (descending)
        by_chr = sorted(
            self.assemblies, key=lambda a: a.stats.chromosome_count, reverse=True
        )

        # Sort by chromosome percentage (descending)
        by_chr_pct = sorted(self.assemblies, key=lambda a: a.chromosome_pct, reverse=True)

        # Sort by fragmentation (ascending - lower is better)
        by_frag = sorted(self.assemblies, key=lambda a: a.fragmentation_index)

        # Sort by total length (descending)
        by_length = sorted(
            self.assemblies, key=lambda a: a.stats.total_length, reverse=True
        )

        return {
            "by_n50": [a.name for a in by_n50],
            "by_chromosome_count": [a.name for a in by_chr],
            "by_chromosome_pct": [a.name for a in by_chr_pct],
            "by_fragmentation": [a.name for a in by_frag],
            "by_total_length": [a.name for a in by_length],
        }


def _check_qc_issues(stats: AssemblyStats, scaffolds: list[ScaffoldInfo]) -> tuple[list[str], list[str]]:
    """Check for QC issues in an assembly.

    Args:
        stats: Assembly statistics
        scaffolds: Classified scaffolds

    Returns:
        Tuple of (warnings, flags)
    """
    warnings = []
    flags = []

    # Check for low chromosome count
    if stats.chromosome_count == 0:
        flags.append("No chromosomes detected")
    elif stats.chromosome_count < 5:
        warnings.append(f"Low chromosome count ({stats.chromosome_count})")

    # Check for high fragmentation
    if stats.total_scaffolds > 0 and stats.chromosome_count > 0:
        frag_ratio = stats.total_scaffolds / stats.chromosome_count
        if frag_ratio > 100:
            flags.append(f"Highly fragmented (ratio: {frag_ratio:.0f})")
        elif frag_ratio > 20:
            warnings.append(f"Moderately fragmented (ratio: {frag_ratio:.0f})")

    # Check for low N50
    if stats.n50 < 1_000_000:
        warnings.append(f"Low N50 ({stats.n50:,} bp)")

    # Check chromosome coverage
    if stats.total_length > 0:
        chr_pct = (stats.chromosome_length / stats.total_length) * 100
        if chr_pct < 50:
            flags.append(f"Low chromosome coverage ({chr_pct:.1f}%)")
        elif chr_pct < 80:
            warnings.append(f"Moderate chromosome coverage ({chr_pct:.1f}%)")

    # Check for unusually large unplaced count
    if stats.unplaced_count > 1000:
        warnings.append(f"Many unplaced scaffolds ({stats.unplaced_count:,})")

    # Check for very small largest scaffold
    if stats.largest_scaffold < 10_000_000:
        warnings.append(f"Small largest scaffold ({stats.largest_scaffold:,} bp)")

    # Check for suspicious chromosome sizes (very uneven)
    chr_scaffolds = [s for s in scaffolds if s.classification == "chromosome"]
    if len(chr_scaffolds) >= 2:
        chr_lengths = [s.length for s in chr_scaffolds]
        size_ratio = max(chr_lengths) / min(chr_lengths)
        if size_ratio > 50:
            warnings.append(f"Uneven chromosome sizes (ratio: {size_ratio:.0f}x)")

    return warnings, flags


def analyze_assembly(
    fasta_path: Path | str,
    min_chromosome_size: int = 10_000_000,
    expected_chromosomes: int | None = None,
) -> AssemblyQC:
    """Analyze a single assembly for QC metrics.

    Args:
        fasta_path: Path to FASTA file
        min_chromosome_size: Minimum size for chromosome classification
        expected_chromosomes: Expected chromosome count

    Returns:
        AssemblyQC with analysis results
    """
    fasta_path = Path(fasta_path)

    # Parse and classify
    scaffolds_raw = parse_fasta(fasta_path)
    scaffolds, stats = classify_scaffolds(
        scaffolds_raw,
        min_chromosome_size=min_chromosome_size,
        expected_chromosomes=expected_chromosomes,
    )

    # Check for QC issues
    warnings, flags = _check_qc_issues(stats, scaffolds)

    return AssemblyQC(
        name=fasta_path.stem,
        path=str(fasta_path),
        stats=stats,
        scaffolds=scaffolds,
        warnings=warnings,
        flags=flags,
    )


def analyze_multiple_assemblies(
    fasta_paths: list[Path | str],
    min_chromosome_size: int = 10_000_000,
    expected_chromosomes: int | None = None,
) -> DashboardResult:
    """Analyze multiple assemblies for comparative QC.

    Args:
        fasta_paths: List of paths to FASTA files
        min_chromosome_size: Minimum size for chromosome classification
        expected_chromosomes: Expected chromosome count (applied to all)

    Returns:
        DashboardResult with all analyses
    """
    assemblies = []

    for path in fasta_paths:
        try:
            qc = analyze_assembly(
                path,
                min_chromosome_size=min_chromosome_size,
                expected_chromosomes=expected_chromosomes,
            )
            assemblies.append(qc)
        except Exception as e:
            # Create a minimal QC entry for failed assemblies
            path = Path(path)
            assemblies.append(
                AssemblyQC(
                    name=path.stem,
                    path=str(path),
                    stats=AssemblyStats(
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
                    ),
                    scaffolds=[],
                    warnings=[],
                    flags=[f"Failed to analyze: {e}"],
                )
            )

    return DashboardResult(assemblies=assemblies)


def format_dashboard_summary(result: DashboardResult) -> str:
    """Format dashboard result as human-readable summary.

    Args:
        result: DashboardResult from analysis

    Returns:
        Formatted summary string
    """
    lines = []
    lines.append("=" * 80)
    lines.append("MULTI-ASSEMBLY QC DASHBOARD")
    lines.append("=" * 80)
    lines.append(f"Generated: {result.generated_at}")
    lines.append(f"Assemblies analyzed: {len(result.assemblies)}")
    lines.append("")

    if result.summary:
        lines.append("SUMMARY STATISTICS:")
        lines.append("-" * 40)
        s = result.summary
        lines.append(f"  N50 range:        {s['n50']['min']:,} - {s['n50']['max']:,} bp")
        lines.append(f"  N50 mean:         {s['n50']['mean']:,} bp")
        lines.append(f"  Best N50:         {s['n50']['best']}")
        lines.append(f"  Chromosome count: {s['chromosome_count']['min']} - {s['chromosome_count']['max']}")
        lines.append(f"  Chr coverage:     {s['chromosome_pct']['min']:.1f}% - {s['chromosome_pct']['max']:.1f}%")
        lines.append(f"  Total warnings:   {s['total_warnings']}")
        lines.append(f"  Total flags:      {s['total_flags']}")
        lines.append("")

    # Rankings
    rankings = result.get_rankings()
    if rankings:
        lines.append("RANKINGS:")
        lines.append("-" * 40)
        lines.append(f"  By N50:           {', '.join(rankings['by_n50'][:5])}")
        lines.append(f"  By chr count:     {', '.join(rankings['by_chromosome_count'][:5])}")
        lines.append(f"  By chr coverage:  {', '.join(rankings['by_chromosome_pct'][:5])}")
        lines.append("")

    # Individual assemblies
    lines.append("ASSEMBLY DETAILS:")
    lines.append("-" * 80)

    for qc in result.assemblies:
        lines.append(f"\n{qc.name}")
        lines.append("  " + "-" * 40)
        lines.append(f"  Scaffolds:        {qc.stats.total_scaffolds:,}")
        lines.append(f"  Total length:     {qc.stats.total_length:,} bp ({qc.stats.total_length/1e9:.2f} Gb)")
        lines.append(f"  N50:              {qc.stats.n50:,} bp ({qc.stats.n50/1e6:.1f} Mb)")
        lines.append(f"  Chromosomes:      {qc.stats.chromosome_count}")
        lines.append(f"  Chr coverage:     {qc.chromosome_pct:.1f}%")
        lines.append(f"  Fragmentation:    {qc.fragmentation_index:.1f}")

        if qc.warnings:
            lines.append(f"  Warnings:         {len(qc.warnings)}")
            for w in qc.warnings:
                lines.append(f"    - {w}")

        if qc.flags:
            lines.append(f"  Flags:            {len(qc.flags)}")
            for f in qc.flags:
                lines.append(f"    ! {f}")

    lines.append("")
    lines.append("=" * 80)

    return "\n".join(lines)


def format_dashboard_tsv(result: DashboardResult) -> str:
    """Format dashboard result as TSV.

    Args:
        result: DashboardResult from analysis

    Returns:
        TSV-formatted string
    """
    lines = []

    # Header
    lines.append("\t".join([
        "name",
        "scaffolds",
        "total_length",
        "n50",
        "chromosomes",
        "chr_length",
        "chr_pct",
        "fragmentation",
        "unplaced",
        "warnings",
        "flags",
    ]))

    # Data rows
    for qc in result.assemblies:
        frag = f"{qc.fragmentation_index:.2f}" if qc.fragmentation_index != float("inf") else "inf"
        lines.append("\t".join([
            qc.name,
            str(qc.stats.total_scaffolds),
            str(qc.stats.total_length),
            str(qc.stats.n50),
            str(qc.stats.chromosome_count),
            str(qc.stats.chromosome_length),
            f"{qc.chromosome_pct:.2f}",
            frag,
            str(qc.stats.unplaced_count),
            str(len(qc.warnings)),
            str(len(qc.flags)),
        ]))

    return "\n".join(lines)


def generate_dashboard_html(result: DashboardResult) -> str:
    """Generate an interactive HTML dashboard.

    Args:
        result: DashboardResult from analysis

    Returns:
        HTML string for the dashboard
    """
    # Prepare data for charts
    names = [a.name for a in result.assemblies]
    n50s = [a.stats.n50 / 1_000_000 for a in result.assemblies]  # Convert to Mb
    chr_counts = [a.stats.chromosome_count for a in result.assemblies]
    chr_pcts = [a.chromosome_pct for a in result.assemblies]
    total_lengths = [a.stats.total_length / 1_000_000_000 for a in result.assemblies]  # Convert to Gb

    # JSON data for JavaScript
    chart_data = json.dumps({
        "names": names,
        "n50s": n50s,
        "chr_counts": chr_counts,
        "chr_pcts": chr_pcts,
        "total_lengths": total_lengths,
    })

    assemblies_json = json.dumps([a.to_dict() for a in result.assemblies])

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChromDetect Multi-Assembly QC Dashboard</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: #f5f5f5;
            color: #333;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        header {{
            background: linear-gradient(135deg, #2c3e50, #3498db);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 20px;
        }}
        header h1 {{
            font-size: 2em;
            margin-bottom: 10px;
        }}
        header p {{
            opacity: 0.9;
        }}
        .summary-cards {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        .card {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .card h3 {{
            font-size: 0.9em;
            color: #666;
            margin-bottom: 5px;
        }}
        .card .value {{
            font-size: 1.8em;
            font-weight: bold;
            color: #2c3e50;
        }}
        .card .detail {{
            font-size: 0.85em;
            color: #888;
        }}
        .charts-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}
        .chart-container {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .chart-container h2 {{
            font-size: 1.1em;
            margin-bottom: 15px;
            color: #2c3e50;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }}
        th {{
            background: #2c3e50;
            color: white;
            font-weight: 600;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .warning {{
            color: #f39c12;
            font-size: 0.9em;
        }}
        .flag {{
            color: #e74c3c;
            font-size: 0.9em;
        }}
        .good {{
            color: #27ae60;
        }}
        .sort-btn {{
            cursor: pointer;
            user-select: none;
        }}
        .sort-btn:hover {{
            text-decoration: underline;
        }}
        footer {{
            text-align: center;
            padding: 20px;
            color: #888;
            font-size: 0.9em;
        }}
        @media (max-width: 768px) {{
            .charts-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Multi-Assembly QC Dashboard</h1>
            <p>Generated: {result.generated_at} | Assemblies: {len(result.assemblies)}</p>
        </header>

        <div class="summary-cards">
            <div class="card">
                <h3>Total Assemblies</h3>
                <div class="value">{len(result.assemblies)}</div>
            </div>
            <div class="card">
                <h3>Best N50</h3>
                <div class="value">{result.summary.get('n50', {}).get('max', 0) / 1_000_000:.1f} Mb</div>
                <div class="detail">{result.summary.get('n50', {}).get('best', 'N/A')}</div>
            </div>
            <div class="card">
                <h3>Max Chromosomes</h3>
                <div class="value">{result.summary.get('chromosome_count', {}).get('max', 0)}</div>
            </div>
            <div class="card">
                <h3>Best Coverage</h3>
                <div class="value">{result.summary.get('chromosome_pct', {}).get('max', 0):.1f}%</div>
            </div>
            <div class="card">
                <h3>Warnings</h3>
                <div class="value {'warning' if result.summary.get('total_warnings', 0) > 0 else 'good'}">{result.summary.get('total_warnings', 0)}</div>
            </div>
            <div class="card">
                <h3>Flags</h3>
                <div class="value {'flag' if result.summary.get('total_flags', 0) > 0 else 'good'}">{result.summary.get('total_flags', 0)}</div>
            </div>
        </div>

        <div class="charts-grid">
            <div class="chart-container">
                <h2>N50 Comparison (Mb)</h2>
                <canvas id="n50Chart"></canvas>
            </div>
            <div class="chart-container">
                <h2>Chromosome Count</h2>
                <canvas id="chrChart"></canvas>
            </div>
            <div class="chart-container">
                <h2>Chromosome Coverage (%)</h2>
                <canvas id="coverageChart"></canvas>
            </div>
            <div class="chart-container">
                <h2>Total Assembly Size (Gb)</h2>
                <canvas id="sizeChart"></canvas>
            </div>
        </div>

        <div class="card" style="margin-bottom: 20px;">
            <h2 style="margin-bottom: 15px;">Assembly Comparison Table</h2>
            <div style="overflow-x: auto;">
                <table id="assemblyTable">
                    <thead>
                        <tr>
                            <th class="sort-btn" onclick="sortTable(0)">Assembly</th>
                            <th class="sort-btn" onclick="sortTable(1)">Scaffolds</th>
                            <th class="sort-btn" onclick="sortTable(2)">Total Size</th>
                            <th class="sort-btn" onclick="sortTable(3)">N50</th>
                            <th class="sort-btn" onclick="sortTable(4)">Chromosomes</th>
                            <th class="sort-btn" onclick="sortTable(5)">Chr %</th>
                            <th>Issues</th>
                        </tr>
                    </thead>
                    <tbody id="tableBody">
                    </tbody>
                </table>
            </div>
        </div>

        <footer>
            Generated by ChromDetect | <a href="https://github.com/your-repo/chromdetect">GitHub</a>
        </footer>
    </div>

    <script>
        const chartData = {chart_data};
        const assemblies = {assemblies_json};

        // Color palette
        const colors = [
            'rgba(52, 152, 219, 0.8)',
            'rgba(46, 204, 113, 0.8)',
            'rgba(155, 89, 182, 0.8)',
            'rgba(241, 196, 15, 0.8)',
            'rgba(230, 126, 34, 0.8)',
            'rgba(231, 76, 60, 0.8)',
            'rgba(26, 188, 156, 0.8)',
            'rgba(52, 73, 94, 0.8)',
        ];

        function getColors(count) {{
            const result = [];
            for (let i = 0; i < count; i++) {{
                result.push(colors[i % colors.length]);
            }}
            return result;
        }}

        // N50 Chart
        new Chart(document.getElementById('n50Chart'), {{
            type: 'bar',
            data: {{
                labels: chartData.names,
                datasets: [{{
                    label: 'N50 (Mb)',
                    data: chartData.n50s,
                    backgroundColor: getColors(chartData.names.length),
                }}]
            }},
            options: {{
                responsive: true,
                plugins: {{ legend: {{ display: false }} }},
                scales: {{ y: {{ beginAtZero: true }} }}
            }}
        }});

        // Chromosome Count Chart
        new Chart(document.getElementById('chrChart'), {{
            type: 'bar',
            data: {{
                labels: chartData.names,
                datasets: [{{
                    label: 'Chromosomes',
                    data: chartData.chr_counts,
                    backgroundColor: getColors(chartData.names.length),
                }}]
            }},
            options: {{
                responsive: true,
                plugins: {{ legend: {{ display: false }} }},
                scales: {{ y: {{ beginAtZero: true }} }}
            }}
        }});

        // Coverage Chart
        new Chart(document.getElementById('coverageChart'), {{
            type: 'bar',
            data: {{
                labels: chartData.names,
                datasets: [{{
                    label: 'Coverage %',
                    data: chartData.chr_pcts,
                    backgroundColor: getColors(chartData.names.length),
                }}]
            }},
            options: {{
                responsive: true,
                plugins: {{ legend: {{ display: false }} }},
                scales: {{ y: {{ beginAtZero: true, max: 100 }} }}
            }}
        }});

        // Size Chart
        new Chart(document.getElementById('sizeChart'), {{
            type: 'bar',
            data: {{
                labels: chartData.names,
                datasets: [{{
                    label: 'Size (Gb)',
                    data: chartData.total_lengths,
                    backgroundColor: getColors(chartData.names.length),
                }}]
            }},
            options: {{
                responsive: true,
                plugins: {{ legend: {{ display: false }} }},
                scales: {{ y: {{ beginAtZero: true }} }}
            }}
        }});

        // Populate table
        function populateTable() {{
            const tbody = document.getElementById('tableBody');
            tbody.innerHTML = '';

            assemblies.forEach(a => {{
                const row = document.createElement('tr');

                const issues = [];
                if (a.warnings.length > 0) {{
                    issues.push(`<span class="warning">${{a.warnings.length}} warnings</span>`);
                }}
                if (a.flags.length > 0) {{
                    issues.push(`<span class="flag">${{a.flags.length}} flags</span>`);
                }}

                row.innerHTML = `
                    <td>${{a.name}}</td>
                    <td>${{a.stats.total_scaffolds.toLocaleString()}}</td>
                    <td>${{(a.stats.total_length / 1e9).toFixed(2)}} Gb</td>
                    <td>${{(a.stats.n50 / 1e6).toFixed(1)}} Mb</td>
                    <td>${{a.stats.chromosome_count}}</td>
                    <td>${{a.chromosome_pct.toFixed(1)}}%</td>
                    <td>${{issues.length > 0 ? issues.join(', ') : '<span class="good">None</span>'}}</td>
                `;
                tbody.appendChild(row);
            }});
        }}

        // Sort table
        let sortDir = {{}};
        function sortTable(col) {{
            const table = document.getElementById('assemblyTable');
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr'));

            sortDir[col] = !sortDir[col];

            rows.sort((a, b) => {{
                let aVal = a.cells[col].textContent;
                let bVal = b.cells[col].textContent;

                // Try to parse as number
                const aNum = parseFloat(aVal.replace(/[^0-9.-]/g, ''));
                const bNum = parseFloat(bVal.replace(/[^0-9.-]/g, ''));

                if (!isNaN(aNum) && !isNaN(bNum)) {{
                    return sortDir[col] ? aNum - bNum : bNum - aNum;
                }}

                return sortDir[col] ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
            }});

            rows.forEach(row => tbody.appendChild(row));
        }}

        populateTable();
    </script>
</body>
</html>"""

    return html
