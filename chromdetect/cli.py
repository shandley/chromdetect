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
from chromdetect.dashboard import (
    analyze_multiple_assemblies,
    format_dashboard_summary,
    format_dashboard_tsv,
    generate_dashboard_html,
)
from chromdetect.html_report import generate_html_report
from chromdetect.karyotype import (
    KaryotypeDatabase,
    format_karyotype_summary,
    format_karyotype_tsv,
    validate_karyotype,
)
from chromdetect.patterns import (
    CHROMOSOME_PATTERNS,
    FRAGMENT_PATTERNS,
    UNLOCALIZED_PATTERNS,
    load_custom_patterns,
    merge_patterns,
)
from chromdetect.standardize import (
    check_ncbi_compliance,
    detect_convention,
    format_compliance_summary,
    format_rename_summary,
    standardize_fasta,
)
from chromdetect.validation import (
    format_validation_summary,
    format_validation_tsv,
    validate_fasta_against_report,
)
from chromdetect.version import (
    compare_fasta_files,
    format_version_json,
    format_version_summary,
    format_version_tsv,
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

Validation mode (check FASTA matches assembly report):
  chromdetect assembly.fasta --assembly-report report.txt --validate
  chromdetect assembly.fasta --assembly-report report.txt --validate --format json
  chromdetect assembly.fasta --assembly-report report.txt --validate --strict

Karyotype validation (check chromosome count matches expected):
  chromdetect assembly.fasta --check-karyotype human
  chromdetect assembly.fasta --check-karyotype "Mus musculus"
  chromdetect --list-species  # List available species

Standardization (rename scaffolds to target convention):
  chromdetect assembly.fasta --rename ucsc -o renamed.fasta
  chromdetect assembly.fasta --rename ensembl -o renamed.fasta
  chromdetect assembly.fasta --rename refseq --species human -o renamed.fasta
  chromdetect assembly.fasta --detect-convention
  chromdetect assembly.fasta --check-compliance

Version comparison (track changes between assembly versions):
  chromdetect v1.fasta --compare-versions v2.fasta
  chromdetect v1.fasta --compare-versions v2.fasta --format json
  chromdetect v1.fasta --compare-versions v2.fasta --format tsv -o changes.tsv

Multi-assembly QC dashboard:
  chromdetect --dashboard *.fasta -o dashboard.html --format html
  chromdetect --dashboard assembly1.fa assembly2.fa assembly3.fa
  chromdetect --dashboard *.fasta --format json -o qc_report.json

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
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate FASTA against assembly report (requires --assembly-report)",
    )
    parser.add_argument(
        "--size-tolerance",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help="Allowed fractional size difference for validation (default: 0.0 = exact)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat validation warnings as errors",
    )
    parser.add_argument(
        "--check-karyotype",
        type=str,
        metavar="SPECIES",
        help="Validate chromosome count/names against expected karyotype for species",
    )
    parser.add_argument(
        "--karyotype-file",
        type=Path,
        metavar="FILE",
        help="Custom karyotype database file (JSON/YAML) for --check-karyotype",
    )
    parser.add_argument(
        "--list-species",
        action="store_true",
        help="List available species in karyotype database and exit",
    )
    parser.add_argument(
        "--rename",
        type=str,
        metavar="CONVENTION",
        choices=["ucsc", "ensembl", "simple", "refseq", "genbank"],
        help="Rename scaffolds to target convention (ucsc, ensembl, simple, refseq, genbank)",
    )
    parser.add_argument(
        "--species",
        type=str,
        default="human",
        metavar="SPECIES",
        help="Species for accession-based renaming (default: human)",
    )
    parser.add_argument(
        "--check-compliance",
        action="store_true",
        help="Check scaffold names for NCBI submission compliance",
    )
    parser.add_argument(
        "--detect-convention",
        action="store_true",
        help="Detect the naming convention used in the FASTA file",
    )
    parser.add_argument(
        "--compare-versions",
        type=Path,
        metavar="FASTA2",
        help="Compare two assembly versions and track scaffold changes",
    )
    parser.add_argument(
        "--version-size-tolerance",
        type=float,
        default=0.1,
        metavar="FLOAT",
        help="Max relative size difference for matching renamed scaffolds (default: 0.1)",
    )
    parser.add_argument(
        "--version-size-threshold",
        type=float,
        default=0.2,
        metavar="FLOAT",
        help="Min relative size change to flag as resized (default: 0.2)",
    )
    parser.add_argument(
        "--dashboard",
        nargs="+",
        type=Path,
        metavar="FASTA",
        help="Generate QC dashboard for multiple assemblies",
    )

    args = parser.parse_args()

    # Handle --list-patterns before requiring fasta argument
    if args.list_patterns:
        show_patterns()
        sys.exit(EXIT_SUCCESS)

    # Handle --list-species before requiring fasta argument
    if args.list_species:
        db = KaryotypeDatabase()
        if args.karyotype_file:
            if not args.karyotype_file.exists():
                print(f"Error: Karyotype file not found: {args.karyotype_file}", file=sys.stderr)
                sys.exit(EXIT_NOINPUT)
            try:
                count = db.load_custom(args.karyotype_file)
                print(f"Loaded {count} custom karyotypes from {args.karyotype_file}", file=sys.stderr)
            except ValueError as e:
                print(f"Error: Invalid karyotype file: {e}", file=sys.stderr)
                sys.exit(EXIT_DATAERR)

        print("Available species in karyotype database:")
        print("-" * 50)
        for species in db.list_species():
            karyotype = db.lookup(species)
            if karyotype:
                print(f"  {species:<20} {karyotype.scientific_name:<30} n={karyotype.chromosome_count}")
        print()
        print(f"Total: {len(db)} species")
        print()
        print("Use --check-karyotype <species> to validate against expected karyotype")
        sys.exit(EXIT_SUCCESS)

    # Handle dashboard mode
    if args.dashboard:
        # Validate all input files exist
        missing = [p for p in args.dashboard if not p.exists()]
        if missing:
            for p in missing:
                print(f"Error: File not found: {p}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        if not args.quiet:
            print(f"Generating QC dashboard for {len(args.dashboard)} assemblies...", file=sys.stderr)

        try:
            dashboard_result = analyze_multiple_assemblies(
                args.dashboard,
                min_chromosome_size=args.min_size,
                expected_chromosomes=args.karyotype,
            )
        except Exception as e:
            print(f"Error during dashboard generation: {e}", file=sys.stderr)
            sys.exit(EXIT_ERROR)

        # Format output
        if args.format == "json":
            output = json.dumps(dashboard_result.to_dict(), indent=2)
        elif args.format == "tsv":
            output = format_dashboard_tsv(dashboard_result)
        elif args.format == "html":
            output = generate_dashboard_html(dashboard_result)
        else:
            output = format_dashboard_summary(dashboard_result)

        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output)
            if not args.quiet:
                print(f"Dashboard written to {args.output}", file=sys.stderr)
        else:
            print(output)

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

    # Handle validation mode
    if args.validate:
        if not args.fasta:
            parser.error("validation requires a FASTA file")
        if not args.assembly_report:
            parser.error("validation requires --assembly-report")

        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"Error: File not found: {args.fasta}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        if not args.assembly_report.exists():
            print(
                f"Error: Assembly report file not found: {args.assembly_report}",
                file=sys.stderr,
            )
            sys.exit(EXIT_NOINPUT)

        if not args.quiet:
            print(f"Validating {args.fasta} against {args.assembly_report}...", file=sys.stderr)

        try:
            result = validate_fasta_against_report(
                fasta_path,
                args.assembly_report,
                size_tolerance=args.size_tolerance,
                strict=args.strict,
            )
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(EXIT_DATAERR)
        except Exception as e:
            print(f"Error during validation: {e}", file=sys.stderr)
            sys.exit(EXIT_ERROR)

        # Format output
        if args.format == "json":
            output = json.dumps(result.to_dict(), indent=2)
        elif args.format == "tsv":
            output = format_validation_tsv(result)
        else:
            output = format_validation_summary(result)

        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output)
            if not args.quiet:
                print(f"Validation results written to {args.output}", file=sys.stderr)
        else:
            print(output)

        # Exit with appropriate code
        sys.exit(EXIT_SUCCESS if result.is_valid else EXIT_ERROR)

    # Handle karyotype validation mode
    if args.check_karyotype:
        if not args.fasta:
            parser.error("karyotype validation requires a FASTA file")

        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"Error: File not found: {args.fasta}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        # Load karyotype database
        karyotype_db = KaryotypeDatabase()
        if args.karyotype_file:
            if not args.karyotype_file.exists():
                print(f"Error: Karyotype file not found: {args.karyotype_file}", file=sys.stderr)
                sys.exit(EXIT_NOINPUT)
            try:
                count = karyotype_db.load_custom(args.karyotype_file)
                if not args.quiet:
                    print(f"Loaded {count} custom karyotypes", file=sys.stderr)
            except ValueError as e:
                print(f"Error: Invalid karyotype file: {e}", file=sys.stderr)
                sys.exit(EXIT_DATAERR)

        if not args.quiet:
            print(f"Parsing {args.fasta}...", file=sys.stderr)

        try:
            scaffolds = parse_fasta(fasta_path)
        except ValueError as e:
            print(f"Error: Invalid FASTA format: {e}", file=sys.stderr)
            sys.exit(EXIT_DATAERR)

        if not args.quiet:
            print(f"Found {len(scaffolds)} scaffolds", file=sys.stderr)
            print("Classifying scaffolds...", file=sys.stderr)

        # Classify scaffolds first to identify chromosomes
        results, stats = classify_scaffolds(
            scaffolds,
            min_chromosome_size=args.min_size,
            expected_chromosomes=args.karyotype,
        )

        # Extract chromosome names/IDs from classification
        chromosome_ids = []
        for r in results:
            if r.classification == "chromosome":
                # Use chromosome_id if available, otherwise use name
                if r.chromosome_id:
                    chromosome_ids.append(r.chromosome_id)
                else:
                    chromosome_ids.append(r.name)

        if not args.quiet:
            print(f"Found {len(chromosome_ids)} chromosome-level scaffolds", file=sys.stderr)
            print(f"Validating against karyotype for '{args.check_karyotype}'...", file=sys.stderr)

        # Validate against karyotype
        karyotype_result = validate_karyotype(
            chromosome_ids,
            args.check_karyotype,
            database=karyotype_db,
            strict=args.strict,
        )

        # Format output
        if args.format == "json":
            output = json.dumps(karyotype_result.to_dict(), indent=2)
        elif args.format == "tsv":
            output = format_karyotype_tsv(karyotype_result)
        else:
            output = format_karyotype_summary(karyotype_result)

        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output)
            if not args.quiet:
                print(f"Karyotype validation results written to {args.output}", file=sys.stderr)
        else:
            print(output)

        sys.exit(EXIT_SUCCESS if karyotype_result.is_valid else EXIT_ERROR)

    # Handle NCBI compliance check
    if args.check_compliance:
        if not args.fasta:
            parser.error("compliance check requires a FASTA file")

        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"Error: File not found: {args.fasta}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        if not args.quiet:
            print(f"Checking NCBI compliance for {args.fasta}...", file=sys.stderr)

        try:
            scaffolds = parse_fasta(fasta_path)
        except ValueError as e:
            print(f"Error: Invalid FASTA format: {e}", file=sys.stderr)
            sys.exit(EXIT_DATAERR)

        compliance_result = check_ncbi_compliance(scaffolds)

        if args.format == "json":
            output = json.dumps(compliance_result.to_dict(), indent=2)
        else:
            output = format_compliance_summary(compliance_result)

        if args.output:
            with open(args.output, "w") as f:
                f.write(output)
            if not args.quiet:
                print(f"Compliance results written to {args.output}", file=sys.stderr)
        else:
            print(output)

        sys.exit(EXIT_SUCCESS if compliance_result.is_compliant else EXIT_ERROR)

    # Handle detect convention
    if args.detect_convention:
        if not args.fasta:
            parser.error("convention detection requires a FASTA file")

        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"Error: File not found: {args.fasta}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        try:
            scaffolds = parse_fasta(fasta_path)
        except ValueError as e:
            print(f"Error: Invalid FASTA format: {e}", file=sys.stderr)
            sys.exit(EXIT_DATAERR)

        names = [s[0] for s in scaffolds]
        detected = detect_convention(names)

        if args.format == "json":
            output = json.dumps({
                "detected_convention": detected.value if detected else None,
                "total_scaffolds": len(names),
                "sample_names": names[:10],
            }, indent=2)
        else:
            if detected:
                output = f"Detected naming convention: {detected.value}\n"
                output += f"Sample names: {', '.join(names[:5])}"
                if len(names) > 5:
                    output += f" ... ({len(names)} total)"
            else:
                output = "Could not determine naming convention\n"
                output += f"Sample names: {', '.join(names[:5])}"

        print(output)
        sys.exit(EXIT_SUCCESS)

    # Handle version comparison
    if args.compare_versions:
        if not args.fasta:
            parser.error("version comparison requires a FASTA file")

        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"Error: File not found: {args.fasta}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        if not args.compare_versions.exists():
            print(f"Error: File not found: {args.compare_versions}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        if not args.quiet:
            print("Comparing assembly versions:", file=sys.stderr)
            print(f"  Version 1: {args.fasta}", file=sys.stderr)
            print(f"  Version 2: {args.compare_versions}", file=sys.stderr)

        try:
            version_result = compare_fasta_files(
                fasta_path,
                args.compare_versions,
                min_chromosome_size=args.min_size,
                size_tolerance=args.version_size_tolerance,
                size_threshold=args.version_size_threshold,
            )
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(EXIT_DATAERR)
        except Exception as e:
            print(f"Error during version comparison: {e}", file=sys.stderr)
            sys.exit(EXIT_ERROR)

        # Format output
        if args.format == "json":
            import json as json_module
            output = json_module.dumps(format_version_json(version_result), indent=2)
        elif args.format == "tsv":
            output = format_version_tsv(version_result)
        else:
            output = format_version_summary(version_result)

        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output)
            if not args.quiet:
                print(f"Version comparison results written to {args.output}", file=sys.stderr)
        else:
            print(output)

        sys.exit(EXIT_SUCCESS)

    # Handle rename/standardize
    if args.rename:
        if not args.fasta:
            parser.error("rename requires a FASTA file")
        if not args.output:
            parser.error("rename requires --output to specify output file")

        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"Error: File not found: {args.fasta}", file=sys.stderr)
            sys.exit(EXIT_NOINPUT)

        if not args.quiet:
            print(f"Renaming scaffolds in {args.fasta} to {args.rename} convention...", file=sys.stderr)

        try:
            rename_result = standardize_fasta(
                fasta_path,
                args.output,
                args.rename,
                species=args.species,
            )
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(EXIT_ERROR)

        if not args.quiet:
            print(f"Renamed {rename_result.renamed_count} scaffolds", file=sys.stderr)
            print(f"Output written to {args.output}", file=sys.stderr)

        # Print summary or JSON
        if args.format == "json":
            print(json.dumps(rename_result.to_dict(), indent=2))
        elif args.verbose:
            print(format_rename_summary(rename_result))

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
    fasta_path: Path | None = None if use_stdin else Path(args.fasta)  # type: ignore[no-redef]

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
