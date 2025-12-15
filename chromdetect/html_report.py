"""
HTML report generation for ChromDetect.

This module generates standalone HTML reports with embedded charts
for visualizing assembly classification results.
"""

from __future__ import annotations

import html
from datetime import datetime

from chromdetect.core import AssemblyStats, ScaffoldInfo

# Import version at module level to avoid circular import
# This is set at the end of the file after the main import
__version__ = "0.5.0"


def _format_bp(bp: int) -> str:
    """Format base pairs in human-readable format."""
    if bp >= 1_000_000_000:
        return f"{bp / 1_000_000_000:.2f} Gb"
    elif bp >= 1_000_000:
        return f"{bp / 1_000_000:.2f} Mb"
    elif bp >= 1_000:
        return f"{bp / 1_000:.2f} kb"
    return f"{bp} bp"


def _generate_pie_chart(
    data: list[tuple[str, int, str]],
    title: str,
    width: int = 300,
    height: int = 250,
) -> str:
    """
    Generate an SVG pie chart.

    Args:
        data: List of (label, value, color) tuples
        title: Chart title
        width: Chart width in pixels
        height: Chart height in pixels

    Returns:
        SVG markup string
    """
    total = sum(v for _, v, _ in data)
    if total == 0:
        return f'<div class="chart-placeholder">{title}: No data</div>'

    cx, cy = width // 2, (height - 30) // 2 + 15
    radius = min(cx, cy) - 40

    svg_parts = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<text x="{width // 2}" y="15" text-anchor="middle" '
        f'font-size="14" font-weight="bold">{html.escape(title)}</text>',
    ]

    start_angle: float = 0.0
    for _label, value, color in data:
        if value == 0:
            continue
        percentage = value / total
        end_angle = start_angle + percentage * 360

        # Calculate arc points
        import math

        start_rad = math.radians(start_angle - 90)
        end_rad = math.radians(end_angle - 90)

        x1 = cx + radius * math.cos(start_rad)
        y1 = cy + radius * math.sin(start_rad)
        x2 = cx + radius * math.cos(end_rad)
        y2 = cy + radius * math.sin(end_rad)

        large_arc = 1 if percentage > 0.5 else 0

        if percentage >= 0.999:
            # Full circle
            svg_parts.append(
                f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="{color}" />'
            )
        else:
            path = (
                f"M {cx},{cy} L {x1},{y1} "
                f"A {radius},{radius} 0 {large_arc},1 {x2},{y2} Z"
            )
            svg_parts.append(f'<path d="{path}" fill="{color}" />')

        start_angle = end_angle

    # Legend
    legend_y = height - 15
    legend_x = 10
    for label, value, color in data:
        if value == 0:
            continue
        percentage = value / total * 100
        svg_parts.append(
            f'<rect x="{legend_x}" y="{legend_y - 8}" width="10" height="10" fill="{color}" />'
        )
        svg_parts.append(
            f'<text x="{legend_x + 14}" y="{legend_y}" font-size="10">'
            f"{html.escape(label)} ({value:,}, {percentage:.1f}%)</text>"
        )
        legend_x += 100

    svg_parts.append("</svg>")
    return "\n".join(svg_parts)


def _generate_bar_chart(
    data: list[tuple[str, int, str]],
    title: str,
    width: int = 600,
    height: int = 300,
    max_bars: int = 20,
) -> str:
    """
    Generate an SVG bar chart.

    Args:
        data: List of (label, value, color) tuples
        title: Chart title
        width: Chart width in pixels
        height: Chart height in pixels
        max_bars: Maximum number of bars to show

    Returns:
        SVG markup string
    """
    if not data:
        return f'<div class="chart-placeholder">{title}: No data</div>'

    data = data[:max_bars]
    max_value = max(v for _, v, _ in data) if data else 1

    chart_left = 80
    chart_top = 35
    chart_width = width - chart_left - 20
    chart_height = height - chart_top - 40
    bar_width = max(10, (chart_width - len(data) * 2) // len(data))

    svg_parts = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<text x="{width // 2}" y="20" text-anchor="middle" '
        f'font-size="14" font-weight="bold">{html.escape(title)}</text>',
    ]

    # Y-axis labels
    for i in range(5):
        y = chart_top + chart_height * (4 - i) / 4
        value = max_value * i / 4
        svg_parts.append(
            f'<text x="{chart_left - 5}" y="{y + 4}" text-anchor="end" font-size="9">'
            f"{_format_bp(int(value))}</text>"
        )
        svg_parts.append(
            f'<line x1="{chart_left}" y1="{y}" x2="{width - 20}" y2="{y}" '
            f'stroke="#eee" stroke-width="1" />'
        )

    # Bars
    x = chart_left + 2
    for label, value, color in data:
        bar_height = (value / max_value) * chart_height if max_value > 0 else 0
        bar_y = chart_top + chart_height - bar_height

        svg_parts.append(
            f'<rect x="{x}" y="{bar_y}" width="{bar_width}" height="{bar_height}" '
            f'fill="{color}" />'
        )

        # Rotated label
        label_y = chart_top + chart_height + 5
        svg_parts.append(
            f'<text x="{x + bar_width // 2}" y="{label_y}" font-size="8" '
            f'transform="rotate(45 {x + bar_width // 2},{label_y})">'
            f"{html.escape(label[:15])}</text>"
        )

        x += bar_width + 2

    svg_parts.append("</svg>")
    return "\n".join(svg_parts)


def generate_html_report(
    results: list[ScaffoldInfo],
    stats: AssemblyStats,
    assembly_name: str = "Assembly",
) -> str:
    """
    Generate a complete HTML report for an assembly analysis.

    Args:
        results: List of ScaffoldInfo from classification
        stats: AssemblyStats summary
        assembly_name: Name to display for the assembly

    Returns:
        Complete HTML document string
    """
    # Count classifications
    chr_count = sum(1 for r in results if r.classification == "chromosome")
    unloc_count = sum(1 for r in results if r.classification == "unlocalized")
    unplaced_count = sum(1 for r in results if r.classification == "unplaced")
    other_count = len(results) - chr_count - unloc_count - unplaced_count

    # Generate pie chart data
    classification_data = [
        ("Chromosome", chr_count, "#4CAF50"),
        ("Unlocalized", unloc_count, "#FF9800"),
        ("Unplaced", unplaced_count, "#9E9E9E"),
        ("Other", other_count, "#607D8B"),
    ]

    # Top scaffolds by size for bar chart
    top_scaffolds = sorted(results, key=lambda r: -r.length)[:20]
    scaffold_size_data = [
        (
            r.name[:20],
            r.length,
            "#4CAF50" if r.classification == "chromosome" else "#9E9E9E",
        )
        for r in top_scaffolds
    ]

    # Classification pie chart
    class_pie = _generate_pie_chart(
        classification_data, "Scaffold Classification", 350, 200
    )

    # Size bar chart
    size_bar = _generate_bar_chart(
        scaffold_size_data, "Top 20 Scaffolds by Size", 700, 300
    )

    # Build scaffold table rows
    table_rows = []
    for r in sorted(results, key=lambda x: -x.length)[:100]:
        gc_str = f"{r.gc_content * 100:.1f}%" if r.gc_content else "N/A"

        class_style = {
            "chromosome": "chr",
            "unlocalized": "unloc",
            "unplaced": "unplaced",
        }.get(r.classification, "other")

        table_rows.append(
            f"""
            <tr>
                <td>{html.escape(r.name)}</td>
                <td class="num">{r.length:,}</td>
                <td><span class="class-{class_style}">{r.classification}</span></td>
                <td class="num">{r.confidence:.2f}</td>
                <td>{r.chromosome_id or "-"}</td>
                <td>{gc_str}</td>
            </tr>
            """
        )

    scaffold_table = "\n".join(table_rows)

    # Pre-format GC content for use in template
    gc_content_str = f"{stats.gc_content * 100:.1f}%" if stats.gc_content else "N/A"

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChromDetect Report: {html.escape(assembly_name)}</title>
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
            padding: 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        header {{
            background: linear-gradient(135deg, #1a237e 0%, #283593 100%);
            color: white;
            padding: 30px;
            border-radius: 8px 8px 0 0;
            margin-bottom: 0;
        }}
        header h1 {{
            font-size: 28px;
            margin-bottom: 5px;
        }}
        header .subtitle {{
            opacity: 0.9;
            font-size: 16px;
        }}
        .content {{
            background: white;
            padding: 30px;
            border-radius: 0 0 8px 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .stat-card {{
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
            text-align: center;
        }}
        .stat-card .value {{
            font-size: 28px;
            font-weight: bold;
            color: #1a237e;
        }}
        .stat-card .label {{
            font-size: 12px;
            color: #666;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .quality-badge {{
            display: inline-block;
            padding: 8px 16px;
            border-radius: 20px;
            font-weight: bold;
            font-size: 14px;
        }}
        .charts {{
            display: grid;
            grid-template-columns: 1fr 2fr;
            gap: 30px;
            margin-bottom: 30px;
        }}
        .chart-container {{
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
            display: flex;
            justify-content: center;
            align-items: center;
        }}
        section {{
            margin-bottom: 30px;
        }}
        section h2 {{
            font-size: 20px;
            color: #1a237e;
            margin-bottom: 15px;
            padding-bottom: 10px;
            border-bottom: 2px solid #e0e0e0;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
        }}
        th, td {{
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        th {{
            background: #f5f5f5;
            font-weight: 600;
            color: #333;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .num {{
            text-align: right;
            font-family: 'Monaco', 'Consolas', monospace;
        }}
        .class-chr {{
            color: #4CAF50;
            font-weight: 600;
        }}
        .class-unloc {{
            color: #FF9800;
        }}
        .class-unplaced {{
            color: #9E9E9E;
        }}
        .class-other {{
            color: #607D8B;
        }}
        .badge {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 10px;
            font-size: 10px;
            font-weight: 600;
            text-transform: uppercase;
        }}
        footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 12px;
        }}
        footer a {{
            color: #1a237e;
        }}
        @media (max-width: 768px) {{
            .charts {{
                grid-template-columns: 1fr;
            }}
            .stats-grid {{
                grid-template-columns: repeat(2, 1fr);
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>ChromDetect Assembly Report</h1>
            <div class="subtitle">{html.escape(assembly_name)} | Generated: {timestamp}</div>
        </header>

        <div class="content">
            <section>
                <h2>Assembly Summary</h2>
                <div class="stats-grid">
                    <div class="stat-card">
                        <div class="value">{stats.total_scaffolds:,}</div>
                        <div class="label">Total Scaffolds</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{_format_bp(stats.total_length)}</div>
                        <div class="label">Total Length</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{_format_bp(stats.n50)}</div>
                        <div class="label">N50</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{stats.chromosome_count}</div>
                        <div class="label">Chromosomes</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{gc_content_str}</div>
                        <div class="label">GC Content</div>
                    </div>
                </div>
            </section>

            <section>
                <h2>Visualizations</h2>
                <div class="charts">
                    <div class="chart-container">
                        {class_pie}
                    </div>
                    <div class="chart-container">
                        {size_bar}
                    </div>
                </div>
            </section>

            <section>
                <h2>Scaffold Details (Top 100)</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Name</th>
                            <th class="num">Length (bp)</th>
                            <th>Classification</th>
                            <th class="num">Confidence</th>
                            <th>Chr ID</th>
                            <th>GC Content</th>
                        </tr>
                    </thead>
                    <tbody>
                        {scaffold_table}
                    </tbody>
                </table>
            </section>

            <section>
                <h2>Classification Statistics</h2>
                <div class="stats-grid">
                    <div class="stat-card">
                        <div class="value">{_format_bp(stats.chromosome_length)}</div>
                        <div class="label">Chromosome Length</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{_format_bp(stats.chromosome_n50)}</div>
                        <div class="label">Chromosome N50</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{stats.unlocalized_count}</div>
                        <div class="label">Unlocalized</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{stats.unplaced_count}</div>
                        <div class="label">Unplaced</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{_format_bp(stats.largest_scaffold)}</div>
                        <div class="label">Largest Scaffold</div>
                    </div>
                    <div class="stat-card">
                        <div class="value">{_format_bp(stats.n90)}</div>
                        <div class="label">N90</div>
                    </div>
                </div>
            </section>
        </div>

        <footer>
            Generated by <a href="https://github.com/shandley/chromdetect">ChromDetect</a> v{__version__}
        </footer>
    </div>
</body>
</html>
"""
    return html_content
