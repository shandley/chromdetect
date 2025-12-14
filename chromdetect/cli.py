"""
Command-line interface for ChromDetect.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from chromdetect.core import (
    ScaffoldInfo,
    AssemblyStats,
    classify_scaffolds,
    parse_fasta,
)


def format_output(
    results: list[ScaffoldInfo],
    stats: AssemblyStats,
    output_format: str = "json",
) -> str:
    """Format results for output.

    Args:
        results: List of ScaffoldInfo from classification
        stats: AssemblyStats summary
        output_format: One of "json", "tsv", or "summary"

    Returns:
        Formatted string for output
    """
    if output_format == "json":
        output = {
            "summary": stats.to_dict(),
            "scaffolds": [r.to_dict() for r in results],
        }
        return json.dumps(output, indent=2)

    elif output_format == "tsv":
        lines = ["name\tlength\tclassification\tconfidence\tmethod\tchromosome_id"]
        for r in results:
            lines.append(
                f"{r.name}\t{r.length}\t{r.classification}\t"
                f"{r.confidence}\t{r.detection_method}\t{r.chromosome_id or ''}"
            )
        return "\n".join(lines)

    elif output_format == "summary":
        lines = [
            "=" * 60,
            "CHROMDETECT ASSEMBLY ANALYSIS",
            "=" * 60,
            "",
            f"Total scaffolds:     {stats.total_scaffolds:,}",
            f"Total length:        {stats.total_length:,} bp ({stats.total_length/1e9:.2f} Gb)",
            f"N50:                 {stats.n50:,} bp ({stats.n50/1e6:.1f} Mb)",
            f"N90:                 {stats.n90:,} bp",
            f"Largest scaffold:    {stats.largest_scaffold:,} bp",
            "",
            "Scaffold Classification:",
            f"  Chromosomes:       {stats.chromosome_count:,} "
            f"({stats.chromosome_length/1e9:.2f} Gb)",
            f"  Unlocalized:       {stats.unlocalized_count:,}",
            f"  Unplaced:          {stats.unplaced_count:,}",
            "",
            f"Chromosome N50:      {stats.chromosome_n50:,} bp "
            f"({stats.chromosome_n50/1e6:.1f} Mb)",
        ]

        if stats.gc_content:
            lines.append(f"GC content:          {stats.gc_content*100:.1f}%")

        lines.extend(["", "-" * 60, "Top 20 Scaffolds:", "-" * 60])

        # Show top scaffolds
        top_results = sorted(results, key=lambda r: -r.length)[:20]
        for r in top_results:
            chr_str = f" ({r.chromosome_id})" if r.chromosome_id else ""
            lines.append(
                f"  {r.name:<30} {r.length:>12,} bp  "
                f"{r.classification:<12} {r.confidence:.2f}{chr_str}"
            )

        return "\n".join(lines)

    else:
        raise ValueError(f"Unknown format: {output_format}")


def main() -> None:
    """Main entry point for chromdetect CLI."""
    parser = argparse.ArgumentParser(
        prog="chromdetect",
        description="Detect chromosome-level scaffolds in genome assemblies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  chromdetect assembly.fasta
  chromdetect assembly.fasta.gz --output results.json
  chromdetect assembly.fasta --karyotype 24 --format summary
  chromdetect assembly.fasta --min-size 5000000 --format tsv > scaffolds.tsv

The tool uses multiple detection strategies:
  - Name patterns (Super_scaffold, Chr, etc.)
  - Size heuristics (large scaffolds likely chromosomes)
  - N50-based detection
  - Optional karyotype-informed adjustment

Supported naming conventions:
  - chr1, chromosome_1, Chr_X
  - Super_scaffold_1, Superscaffold_1, SUPER_1
  - LG_1 (linkage groups)
  - NC_*, CM* (NCBI accessions)
  - HiC_scaffold_1, Scaffold_1_RaGOO
        """,
    )

    parser.add_argument(
        "fasta",
        type=Path,
        help="Input FASTA file (can be gzipped)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output file (default: stdout)",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["json", "tsv", "summary"],
        default="summary",
        help="Output format (default: summary)",
    )
    parser.add_argument(
        "-k",
        "--karyotype",
        type=int,
        help="Expected chromosome count (for karyotype-informed detection)",
    )
    parser.add_argument(
        "-s",
        "--min-size",
        type=int,
        default=10_000_000,
        help="Minimum size (bp) to consider chromosome-level (default: 10Mb)",
    )
    parser.add_argument(
        "-c",
        "--chromosomes-only",
        action="store_true",
        help="Only output chromosome-level scaffolds",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress progress messages",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 0.1.0",
    )

    args = parser.parse_args()

    if not args.fasta.exists():
        print(f"Error: File not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"Parsing {args.fasta}...", file=sys.stderr)

    try:
        scaffolds = parse_fasta(args.fasta)
    except Exception as e:
        print(f"Error parsing FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"Found {len(scaffolds)} scaffolds", file=sys.stderr)
        print("Classifying scaffolds...", file=sys.stderr)

    results, stats = classify_scaffolds(
        scaffolds,
        min_chromosome_size=args.min_size,
        expected_chromosomes=args.karyotype,
    )

    # Filter if requested
    if args.chromosomes_only:
        results = [r for r in results if r.classification == "chromosome"]

    # Format output
    output = format_output(results, stats, args.format)

    # Write output
    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        if not args.quiet:
            print(f"Results written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
