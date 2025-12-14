"""
Command-line interface for ChromDetect.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from chromdetect import __version__
from chromdetect.core import (
    AssemblyStats,
    ScaffoldInfo,
    classify_scaffolds,
    parse_fasta,
    parse_fasta_from_handle,
)
from chromdetect.patterns import (
    CHROMOSOME_PATTERNS,
    FRAGMENT_PATTERNS,
    UNLOCALIZED_PATTERNS,
)

# Exit codes following sysexits.h conventions
EXIT_SUCCESS = 0
EXIT_ERROR = 1  # General error


def show_patterns() -> None:
    """Display supported naming patterns."""
    print("ChromDetect Supported Naming Patterns")
    print("=" * 50)
    print()
    print("CHROMOSOME PATTERNS (detected as chromosome-level):")
    print("-" * 50)
    for pattern, name in CHROMOSOME_PATTERNS:
        # Clean up regex for display
        display_pattern = pattern.replace("^", "").replace("$", "")
        print(f"  {name:<25} {display_pattern}")

    print()
    print("UNLOCALIZED PATTERNS (chromosome-associated but not placed):")
    print("-" * 50)
    for pattern in UNLOCALIZED_PATTERNS:
        print(f"  {pattern}")

    print()
    print("FRAGMENT PATTERNS (contigs/fragments, not chromosome-level):")
    print("-" * 50)
    for pattern in FRAGMENT_PATTERNS:
        print(f"  {pattern}")


EXIT_USAGE = 2  # Command line usage error (argparse default)
EXIT_NOINPUT = 66  # Input file not found or not readable
EXIT_DATAERR = 65  # Input data format error


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
        type=str,
        nargs="?",
        help="Input FASTA file (can be gzipped), or '-' for stdin",
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
        "--min-confidence",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help="Minimum confidence threshold (0.0-1.0) to include scaffolds",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=0,
        metavar="BP",
        help="Minimum scaffold length (bp) to include in output",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress progress messages",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show detailed processing information",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--list-patterns",
        action="store_true",
        help="Show supported naming patterns and exit",
    )

    args = parser.parse_args()

    # Handle --list-patterns before requiring fasta argument
    if args.list_patterns:
        show_patterns()
        sys.exit(EXIT_SUCCESS)

    # Require fasta argument for normal operation
    if args.fasta is None:
        parser.error("the following arguments are required: fasta")

    # Handle conflicting verbosity flags
    verbose = args.verbose and not args.quiet

    # Handle stdin vs file input
    use_stdin = args.fasta == "-"
    fasta_path = None if use_stdin else Path(args.fasta)

    if verbose:
        print(f"ChromDetect {__version__}", file=sys.stderr)
        print(f"Input: {'stdin' if use_stdin else args.fasta}", file=sys.stderr)
        print(f"Output format: {args.format}", file=sys.stderr)
        print(f"Min chromosome size: {args.min_size:,} bp", file=sys.stderr)
        if args.karyotype:
            print(f"Expected karyotype: {args.karyotype}", file=sys.stderr)

    if not use_stdin:
        assert fasta_path is not None  # for type checker
        if not fasta_path.exists():
            print(f"Error: File not found: {args.fasta}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

    if not args.quiet:
        source = "stdin" if use_stdin else str(fasta_path)
        print(f"Parsing {source}...", file=sys.stderr)

    try:
        if use_stdin:
            scaffolds = parse_fasta_from_handle(sys.stdin)
        else:
            assert fasta_path is not None  # for type checker
            scaffolds = parse_fasta(fasta_path)
    except ValueError as e:
        print(f"Error: Invalid FASTA format: {e}", file=sys.stderr)
        sys.exit(EXIT_DATAERR)
    except OSError as e:
        print(f"Error: Cannot read file: {e}", file=sys.stderr)
        sys.exit(EXIT_NOINPUT)
    except Exception as e:
        print(f"Error parsing FASTA: {e}", file=sys.stderr)
        sys.exit(EXIT_ERROR)

    if not args.quiet:
        print(f"Found {len(scaffolds)} scaffolds", file=sys.stderr)
        print("Classifying scaffolds...", file=sys.stderr)

    if verbose:
        total_length = sum(s[1] for s in scaffolds)
        print(f"Total assembly length: {total_length:,} bp", file=sys.stderr)
        lengths = sorted([s[1] for s in scaffolds], reverse=True)
        print(f"Largest scaffold: {lengths[0]:,} bp", file=sys.stderr)
        print(f"Smallest scaffold: {lengths[-1]:,} bp", file=sys.stderr)

    results, stats = classify_scaffolds(
        scaffolds,
        min_chromosome_size=args.min_size,
        expected_chromosomes=args.karyotype,
    )

    if verbose:
        print("Classification complete:", file=sys.stderr)
        print(f"  Chromosomes: {stats.chromosome_count}", file=sys.stderr)
        print(f"  Unlocalized: {stats.unlocalized_count}", file=sys.stderr)
        print(f"  Unplaced: {stats.unplaced_count}", file=sys.stderr)
        print(f"  N50: {stats.n50:,} bp", file=sys.stderr)

    # Apply filters
    original_count = len(results)

    if args.chromosomes_only:
        results = [r for r in results if r.classification == "chromosome"]

    if args.min_confidence > 0:
        results = [r for r in results if r.confidence >= args.min_confidence]

    if args.min_length > 0:
        results = [r for r in results if r.length >= args.min_length]

    if verbose and len(results) != original_count:
        filters_applied = []
        if args.chromosomes_only:
            filters_applied.append("chromosomes-only")
        if args.min_confidence > 0:
            filters_applied.append(f"min-confidence={args.min_confidence}")
        if args.min_length > 0:
            filters_applied.append(f"min-length={args.min_length:,}")
        print(
            f"Filtered to {len(results)} scaffolds "
            f"(from {original_count} total, filters: {', '.join(filters_applied)})",
            file=sys.stderr,
        )

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
