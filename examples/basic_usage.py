#!/usr/bin/env python3
"""
Basic usage example for ChromDetect.

This example demonstrates how to use ChromDetect to classify scaffolds
in a genome assembly.
"""

from pathlib import Path

from chromdetect import classify_scaffolds, parse_fasta


def main() -> None:
    """Run basic classification example."""
    # Example with a local FASTA file
    # Replace with your actual assembly path
    assembly_path = Path("assembly.fasta")

    if not assembly_path.exists():
        print("Example: Creating a mock assembly for demonstration")
        create_mock_assembly(assembly_path)

    print(f"Parsing {assembly_path}...")
    scaffolds = parse_fasta(assembly_path)
    print(f"Found {len(scaffolds)} scaffolds")

    # Classify scaffolds
    print("\nClassifying scaffolds...")
    results, stats = classify_scaffolds(scaffolds)

    # Print summary statistics
    print("\n" + "=" * 50)
    print("ASSEMBLY SUMMARY")
    print("=" * 50)
    print(f"Total scaffolds:     {stats.total_scaffolds:,}")
    print(f"Total length:        {stats.total_length:,} bp")
    print(f"N50:                 {stats.n50:,} bp")
    print(f"Largest scaffold:    {stats.largest_scaffold:,} bp")
    print()
    print("Classification breakdown:")
    print(f"  Chromosomes:       {stats.chromosome_count}")
    print(f"  Unlocalized:       {stats.unlocalized_count}")
    print(f"  Unplaced:          {stats.unplaced_count}")

    # Print chromosome details
    print("\n" + "-" * 50)
    print("DETECTED CHROMOSOMES")
    print("-" * 50)

    chromosomes = [r for r in results if r.classification == "chromosome"]
    chromosomes.sort(key=lambda r: -r.length)

    for r in chromosomes[:10]:  # Top 10
        chr_id = f" ({r.chromosome_id})" if r.chromosome_id else ""
        print(
            f"  {r.name:<30} {r.length:>12,} bp  "
            f"conf={r.confidence:.2f}{chr_id}"
        )

    if len(chromosomes) > 10:
        print(f"  ... and {len(chromosomes) - 10} more")

    # Clean up mock file
    if assembly_path.name == "assembly.fasta":
        assembly_path.unlink()


def create_mock_assembly(path: Path) -> None:
    """Create a mock assembly file for demonstration."""
    with open(path, "w") as f:
        # Chromosome-level scaffolds
        for i in range(1, 6):
            f.write(f">chr{i}\n")
            f.write("A" * (50_000_000 - i * 5_000_000) + "\n")

        # Super scaffolds
        for i in range(1, 4):
            f.write(f">Super_scaffold_{i}\n")
            f.write("G" * (20_000_000 - i * 3_000_000) + "\n")

        # Unplaced scaffolds
        for i in range(1, 20):
            f.write(f">scaffold_arrow_ctg{i}\n")
            f.write("C" * (500_000 - i * 20_000) + "\n")


if __name__ == "__main__":
    main()
