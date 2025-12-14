#!/usr/bin/env python3
"""
Karyotype-informed classification example.

This example shows how to use known karyotype information to improve
chromosome detection accuracy.
"""

import tempfile
from pathlib import Path

from chromdetect import classify_scaffolds, parse_fasta


def main() -> None:
    """Demonstrate karyotype-informed classification."""

    # Create a mock assembly with inconsistent naming
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as f:
        # Only 2 scaffolds have chromosome names
        f.write(">chr1\n")
        f.write("A" * 100_000_000 + "\n")
        f.write(">chr2\n")
        f.write("G" * 90_000_000 + "\n")

        # But 3 more are clearly chromosome-level by size
        f.write(">scaffold_001\n")
        f.write("C" * 80_000_000 + "\n")
        f.write(">scaffold_002\n")
        f.write("T" * 70_000_000 + "\n")
        f.write(">scaffold_003\n")
        f.write("A" * 60_000_000 + "\n")

        # Small unplaced scaffolds
        for i in range(10):
            f.write(f">contig_{i:03d}\n")
            f.write("G" * 100_000 + "\n")

        assembly_path = Path(f.name)

    print("=" * 60)
    print("KARYOTYPE-INFORMED CLASSIFICATION DEMO")
    print("=" * 60)

    # Parse assembly
    scaffolds = parse_fasta(assembly_path)

    # Without karyotype information
    print("\n1. Classification WITHOUT karyotype information:")
    print("-" * 50)

    results_no_karyo, stats_no_karyo = classify_scaffolds(scaffolds)
    print(f"   Chromosomes detected: {stats_no_karyo.chromosome_count}")

    print("\n   Classified as chromosome:")
    for r in results_no_karyo:
        if r.classification == "chromosome":
            print(f"     {r.name:<20} {r.length:>12,} bp  ({r.detection_method})")

    # With karyotype = 5 (we know there are 5 chromosomes)
    print("\n2. Classification WITH karyotype=5:")
    print("-" * 50)

    results_karyo, stats_karyo = classify_scaffolds(
        scaffolds, expected_chromosomes=5
    )
    print(f"   Chromosomes detected: {stats_karyo.chromosome_count}")

    print("\n   Classified as chromosome:")
    for r in results_karyo:
        if r.classification == "chromosome":
            print(f"     {r.name:<20} {r.length:>12,} bp  ({r.detection_method})")

    # Comparison
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Without karyotype: {stats_no_karyo.chromosome_count} chromosomes")
    print(f"With karyotype=5:  {stats_karyo.chromosome_count} chromosomes")

    if stats_karyo.chromosome_count > stats_no_karyo.chromosome_count:
        promoted = stats_karyo.chromosome_count - stats_no_karyo.chromosome_count
        print(f"\nKaryotype info promoted {promoted} scaffold(s) to chromosome status")

    # Clean up
    assembly_path.unlink()


if __name__ == "__main__":
    main()
