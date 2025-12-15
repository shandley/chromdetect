"""
Command-line interface for ChromDetect.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from chromdetect import __version__
from chromdetect.assembly_report import parse_assembly_report
from chromdetect.compare import (
    compare_assemblies,
    format_comparison_summary,
    format_comparison_tsv,
)
from chromdetect.core import (
    AssemblyStats,
    ScaffoldInfo,
    classify_scaffolds,
    format_bed,
    format_gff,
    parse_fasta,
    parse_fasta_from_handle,
    write_fasta,
)
from chromdetect.html_report import generate_html_report
from chromdetect.patterns import (
    CHROMOSOME_PATTERNS,
    FRAGMENT_PATTERNS,
    UNLOCALIZED_PATTERNS,
    load_custom_patterns,
    merge_patterns,
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
    assembly_name: str = "Assembly",
) -> str:
    """Format results for output.

    Args:
        results: List of ScaffoldInfo from classification
        stats: AssemblyStats summary
        output_format: One of "json", "tsv", "summary", "bed", "gff", or "html"
        assembly_name: Name of the assembly for HTML reports

    Returns:
        Formatted string for output
    """
    if output_format == "html":
        return generate_html_report(results, stats, assembly_name)

    if output_format == "json":
        output = {
            "summary": stats.to_dict(),
            "scaffolds": [r.to_dict() for r in results],
        }
        return json.dumps(output, indent=2)

    elif output_format == "tsv":
        lines = ["name\tlength\tclassification\tconfidence\tmethod\tchromosome_id\tgc_content"]
        for r in results:
            gc_str = f"{r.gc_content:.4f}" if r.gc_content is not None else ""
            lines.append(
                f"{r.name}\t{r.length}\t{r.classification}\t"
                f"{r.confidence}\t{r.detection_method}\t{r.chromosome_id or ''}\t{gc_str}"
            )
        return "\n".join(lines)

    elif output_format == "bed":
        return format_bed(results)

    elif output_format == "gff":
        return format_gff(results)

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
            gc_str = f" GC:{r.gc_content*100:.1f}%" if r.gc_content is not None else ""
            lines.append(
                f"  {r.name:<30} {r.length:>12,} bp  "
                f"{r.classification:<12} {r.confidence:.2f}{chr_str}{gc_str}"
            )

        return "\n".join(lines)

    else:
        raise ValueError(f"Unknown format: {output_format}")


def process_batch(
    args: argparse.Namespace,
    custom_patterns: tuple | None = None,
    assembly_report: object | None = None,
) -> None:
    """Process all FASTA files in a directory.

    Args:
        args: Parsed command-line arguments
        custom_patterns: Optional custom patterns tuple from merge_patterns()
        assembly_report: Optional AssemblyReport for scaffold classification
    """
    batch_dir = args.batch
    if not batch_dir.is_dir():
        print(f"Error: {batch_dir} is not a directory", file=sys.stderr)
        sys.exit(EXIT_NOINPUT)

    # Find all FASTA files
    fasta_extensions = {".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz"}
    fasta_files = []
    for ext in fasta_extensions:
        if ext.endswith(".gz"):
            # Handle double extensions like .fasta.gz
            base_ext = ext[:-3]  # Remove .gz
            fasta_files.extend(batch_dir.glob(f"*{base_ext}.gz"))
        else:
            fasta_files.extend(batch_dir.glob(f"*{ext}"))

    # Remove duplicates and sort
    fasta_files = sorted(set(fasta_files))

    if not fasta_files:
        print(f"Error: No FASTA files found in {batch_dir}", file=sys.stderr)
        sys.exit(EXIT_NOINPUT)

    # Determine output directory
    output_dir = args.output if args.output else batch_dir / "chromdetect_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    verbose = args.verbose and not args.quiet

    if not args.quiet:
        print(f"Found {len(fasta_files)} FASTA files in {batch_dir}", file=sys.stderr)
        print(f"Results will be written to {output_dir}", file=sys.stderr)

    # Determine output file extension based on format
    format_extensions = {
        "json": ".json",
        "tsv": ".tsv",
        "summary": ".txt",
        "bed": ".bed",
        "gff": ".gff",
    }
    out_ext = format_extensions.get(args.format, ".txt")

    # Process each file
    results_summary = []
    for i, fasta_path in enumerate(fasta_files, 1):
        if not args.quiet:
            print(f"[{i}/{len(fasta_files)}] Processing {fasta_path.name}...", file=sys.stderr)

        try:
            # Parse FASTA
            need_full_sequences = args.extract_chromosomes is not None
            scaffolds = parse_fasta(fasta_path, need_full_sequences)

            # Classify
            results, stats = classify_scaffolds(
                scaffolds,
                min_chromosome_size=args.min_size,
                expected_chromosomes=args.karyotype,
                custom_patterns=custom_patterns,
                assembly_report=assembly_report,
            )

            # Handle chromosome extraction if requested
            if args.extract_chromosomes:
                scaffold_seqs = {name: seq for name, _length, seq in scaffolds}
                chr_names = {r.name for r in results if r.classification == "chromosome"}
                chr_sequences = [
                    (name, scaffold_seqs[name]) for name in chr_names if name in scaffold_seqs
                ]
                chr_sequences.sort(key=lambda x: x[0])
                if chr_sequences:
                    # Create output filename based on input
                    chr_out = output_dir / f"{fasta_path.stem}_chromosomes.fasta"
                    write_fasta(chr_sequences, chr_out)

            # Apply filters
            if args.chromosomes_only:
                results = [r for r in results if r.classification == "chromosome"]
            if args.min_confidence > 0:
                results = [r for r in results if r.confidence >= args.min_confidence]
            if args.min_length > 0:
                results = [r for r in results if r.length >= args.min_length]

            # Format and write output
            output = format_output(results, stats, args.format, fasta_path.stem)
            out_ext_actual = ".html" if args.format == "html" else out_ext
            out_file = output_dir / f"{fasta_path.stem}{out_ext_actual}"
            with open(out_file, "w") as f:
                f.write(output)

            # Track summary
            results_summary.append({
                "file": fasta_path.name,
                "scaffolds": stats.total_scaffolds,
                "chromosomes": stats.chromosome_count,
                "total_length": stats.total_length,
                "n50": stats.n50,
            })

            if verbose:
                print(
                    f"  -> {stats.chromosome_count} chromosomes, "
                    f"N50={stats.n50:,} bp",
                    file=sys.stderr,
                )

        except Exception as e:
            print(f"  Error processing {fasta_path.name}: {e}", file=sys.stderr)
            results_summary.append({
                "file": fasta_path.name,
                "error": str(e),
            })

    # Write batch summary
    summary_file = output_dir / "batch_summary.tsv"
    with open(summary_file, "w") as f:
        f.write("file\tscaffolds\tchromosomes\ttotal_length\tn50\terror\n")
        for item in results_summary:
            if "error" in item:
                f.write(f"{item['file']}\t\t\t\t\t{item['error']}\n")
            else:
                f.write(
                    f"{item['file']}\t{item['scaffolds']}\t{item['chromosomes']}\t"
                    f"{item['total_length']}\t{item['n50']}\t\n"
                )

    if not args.quiet:
        successful = sum(1 for r in results_summary if "error" not in r)
        print(f"\nBatch complete: {successful}/{len(fasta_files)} files processed", file=sys.stderr)
        print(f"Summary written to {summary_file}", file=sys.stderr)


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
  chromdetect assembly.fasta --format bed > scaffolds.bed
  chromdetect assembly.fasta --format gff > scaffolds.gff
  chromdetect assembly.fasta --extract-chromosomes -o chromosomes.fasta
  chromdetect --batch assemblies/ --output results/

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
        choices=["json", "tsv", "summary", "bed", "gff", "html"],
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
        "--extract-chromosomes",
        type=Path,
        metavar="FILE",
        help="Extract chromosome sequences to FASTA file",
    )
    parser.add_argument(
        "--batch",
        type=Path,
        metavar="DIR",
        help="Process all FASTA files in directory (batch mode)",
    )
    parser.add_argument(
        "--compare",
        type=Path,
        metavar="FASTA2",
        help="Compare with a second assembly (side-by-side analysis)",
    )
    parser.add_argument(
        "--patterns",
        type=Path,
        metavar="FILE",
        help="Custom patterns file (YAML or JSON) for scaffold name matching",
    )
    parser.add_argument(
        "--assembly-report",
        type=Path,
        metavar="FILE",
        help="NCBI assembly report file for authoritative scaffold classification",
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

    # Handle batch mode
    if args.batch:
        # Load custom patterns for batch mode
        batch_patterns = None
        if args.patterns:
            if not args.patterns.exists():
                print(f"Error: Patterns file not found: {args.patterns}", file=sys.stderr)
                sys.exit(EXIT_NOINPUT)
            try:
                chr_pats, unloc_pats, frag_pats = load_custom_patterns(args.patterns)
                batch_patterns = merge_patterns(chr_pats, unloc_pats, frag_pats)
            except ValueError as e:
                print(f"Error: Invalid patterns file: {e}", file=sys.stderr)
                sys.exit(EXIT_DATAERR)
        # Load assembly report for batch mode
        batch_report = None
        if args.assembly_report:
            if not args.assembly_report.exists():
                print(
                    f"Error: Assembly report file not found: {args.assembly_report}",
                    file=sys.stderr,
                )
                sys.exit(EXIT_NOINPUT)
            try:
                batch_report = parse_assembly_report(args.assembly_report)
            except ValueError as e:
                print(f"Error: Invalid assembly report: {e}", file=sys.stderr)
                sys.exit(EXIT_DATAERR)
        process_batch(args, batch_patterns, batch_report)
        sys.exit(EXIT_SUCCESS)

    # Require fasta argument for normal operation
    if args.fasta is None:
        parser.error("the following arguments are required: fasta")

    # Handle conflicting verbosity flags
    verbose = args.verbose and not args.quiet

    # Load custom patterns if provided
    custom_patterns = None
    if args.patterns:
        if not args.patterns.exists():
            print(f"Error: Patterns file not found: {args.patterns}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)
        try:
            chr_pats, unloc_pats, frag_pats = load_custom_patterns(args.patterns)
            custom_patterns = merge_patterns(chr_pats, unloc_pats, frag_pats)
            if verbose:
                print(
                    f"Loaded {len(chr_pats)} custom chromosome patterns, "
                    f"{len(unloc_pats)} unlocalized patterns, "
                    f"{len(frag_pats)} fragment patterns",
                    file=sys.stderr,
                )
        except ValueError as e:
            print(f"Error: Invalid patterns file: {e}", file=sys.stderr)
            sys.exit(EXIT_DATAERR)

    # Load assembly report if provided
    assembly_report = None
    if args.assembly_report:
        if not args.assembly_report.exists():
            print(
                f"Error: Assembly report file not found: {args.assembly_report}",
                file=sys.stderr,
            )
            sys.exit(EXIT_NOINPUT)
        try:
            assembly_report = parse_assembly_report(args.assembly_report)
            if verbose:
                print(
                    f"Loaded assembly report: {assembly_report.assembly_name or 'unknown'}",
                    file=sys.stderr,
                )
                print(
                    f"  {len(assembly_report.entries)} sequences, "
                    f"{assembly_report.get_expected_chromosome_count()} chromosomes",
                    file=sys.stderr,
                )
        except ValueError as e:
            print(f"Error: Invalid assembly report: {e}", file=sys.stderr)
            sys.exit(EXIT_DATAERR)

    # Handle stdin vs file input
    use_stdin = args.fasta == "-"
    fasta_path = None if use_stdin else Path(args.fasta)

    # Determine if we need full sequences (for extraction)
    need_full_sequences = args.extract_chromosomes is not None

    if verbose:
        print(f"ChromDetect {__version__}", file=sys.stderr)
        print(f"Input: {'stdin' if use_stdin else args.fasta}", file=sys.stderr)
        print(f"Output format: {args.format}", file=sys.stderr)
        print(f"Min chromosome size: {args.min_size:,} bp", file=sys.stderr)
        if args.karyotype:
            print(f"Expected karyotype: {args.karyotype}", file=sys.stderr)
        if args.extract_chromosomes:
            print(f"Extract chromosomes to: {args.extract_chromosomes}", file=sys.stderr)

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
            scaffolds = parse_fasta_from_handle(sys.stdin, need_full_sequences)
        else:
            assert fasta_path is not None  # for type checker
            scaffolds = parse_fasta(fasta_path, need_full_sequences)
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
        custom_patterns=custom_patterns,
        assembly_report=assembly_report,
    )

    if verbose:
        print("Classification complete:", file=sys.stderr)
        print(f"  Chromosomes: {stats.chromosome_count}", file=sys.stderr)
        print(f"  Unlocalized: {stats.unlocalized_count}", file=sys.stderr)
        print(f"  Unplaced: {stats.unplaced_count}", file=sys.stderr)
        print(f"  N50: {stats.n50:,} bp", file=sys.stderr)

    # Handle chromosome extraction
    if args.extract_chromosomes:
        # Build a mapping of scaffold name to sequence
        scaffold_seqs = {name: seq for name, _length, seq in scaffolds}
        # Get chromosome names
        chr_names = {r.name for r in results if r.classification == "chromosome"}
        # Extract chromosome sequences
        chr_sequences = [
            (name, scaffold_seqs[name]) for name in chr_names if name in scaffold_seqs
        ]
        # Sort by name for consistent output
        chr_sequences.sort(key=lambda x: x[0])

        if not chr_sequences:
            print("Warning: No chromosome sequences to extract", file=sys.stderr)
        else:
            write_fasta(chr_sequences, args.extract_chromosomes)
            if not args.quiet:
                print(
                    f"Extracted {len(chr_sequences)} chromosome sequences "
                    f"to {args.extract_chromosomes}",
                    file=sys.stderr,
                )

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

    # Handle comparison mode
    if args.compare:
        if not args.compare.exists():
            print(f"Error: Comparison file not found: {args.compare}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        if not args.quiet:
            print(f"Parsing second assembly {args.compare}...", file=sys.stderr)

        try:
            scaffolds2 = parse_fasta(args.compare, need_full_sequences)
        except ValueError as e:
            print(f"Error: Invalid FASTA format in comparison file: {e}", file=sys.stderr)
            sys.exit(EXIT_DATAERR)
        except OSError as e:
            print(f"Error: Cannot read comparison file: {e}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        if not args.quiet:
            print(f"Found {len(scaffolds2)} scaffolds in second assembly", file=sys.stderr)
            print("Classifying second assembly...", file=sys.stderr)

        results2, stats2 = classify_scaffolds(
            scaffolds2,
            min_chromosome_size=args.min_size,
            expected_chromosomes=args.karyotype,
            custom_patterns=custom_patterns,
            assembly_report=assembly_report,
        )

        # Get assembly names from file paths
        assembly1_name = Path(args.fasta).stem if args.fasta != "-" else "Assembly 1"
        assembly2_name = args.compare.stem

        comparison = compare_assemblies(
            results, stats, results2, stats2,
            assembly1_name=assembly1_name,
            assembly2_name=assembly2_name,
        )

        # Format comparison output
        if args.format == "json":
            output = json.dumps(comparison.to_dict(), indent=2)
        elif args.format == "tsv":
            output = format_comparison_tsv(comparison)
        else:
            output = format_comparison_summary(comparison)

        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output)
            if not args.quiet:
                print(f"Comparison results written to {args.output}", file=sys.stderr)
        else:
            print(output)

        sys.exit(EXIT_SUCCESS)

    # Format output
    assembly_name = Path(args.fasta).stem if args.fasta != "-" else "Assembly"
    output = format_output(results, stats, args.format, assembly_name)

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
