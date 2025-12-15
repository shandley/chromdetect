#!/usr/bin/env python3
"""
Batch processing example for ChromDetect.

This example shows how to process multiple assemblies and compare results.
"""

import json
from pathlib import Path

from chromdetect import classify_scaffolds, parse_fasta

# Example species with known karyotypes
SPECIES_INFO = {
    "human": {"file": "GRCh38.fasta", "karyotype": 23},
    "mouse": {"file": "GRCm39.fasta", "karyotype": 20},
    "zebrafish": {"file": "GRCz11.fasta", "karyotype": 25},
    "fruit_fly": {"file": "dm6.fasta", "karyotype": 4},
}


def process_assembly(fasta_path: Path, karyotype: int | None = None) -> dict:
    """
    Process a single assembly and return results.

    Args:
        fasta_path: Path to FASTA file
        karyotype: Expected chromosome count (optional)

    Returns:
        Dictionary with classification results
    """
    scaffolds = parse_fasta(fasta_path)
    results, stats = classify_scaffolds(
        scaffolds,
        expected_chromosomes=karyotype,
    )

    return {
        "file": str(fasta_path),
        "total_scaffolds": stats.total_scaffolds,
        "total_length_gb": round(stats.total_length / 1e9, 2),
        "n50_mb": round(stats.n50 / 1e6, 1),
        "chromosomes_detected": stats.chromosome_count,
        "chromosomes_expected": karyotype,
        "match": stats.chromosome_count == karyotype if karyotype else None,
        "chromosome_coverage_pct": round(
            100 * stats.chromosome_length / stats.total_length, 1
        ),
    }


def batch_process(assembly_dir: Path) -> list[dict]:
    """
    Process all assemblies in a directory.

    Args:
        assembly_dir: Directory containing FASTA files

    Returns:
        List of result dictionaries
    """
    results = []

    for name, info in SPECIES_INFO.items():
        fasta_path = assembly_dir / info["file"]

        if fasta_path.exists():
            print(f"Processing {name}...")
            result = process_assembly(fasta_path, info["karyotype"])
            result["species"] = name
            results.append(result)
        else:
            print(f"Skipping {name} - file not found: {fasta_path}")

    return results


def print_summary(results: list[dict]) -> None:
    """Print a summary table of results."""
    print("\n" + "=" * 80)
    print("BATCH PROCESSING RESULTS")
    print("=" * 80)

    # Header
    print(
        f"{'Species':<15} {'Size (Gb)':<10} {'N50 (Mb)':<10} "
        f"{'Chr Found':<10} {'Expected':<10} {'Match':<8}"
    )
    print("-" * 80)

    # Data rows
    for r in results:
        match_str = "Yes" if r.get("match") else "No" if r.get("match") is False else "-"
        expected = r.get("chromosomes_expected", "-")
        print(
            f"{r.get('species', 'unknown'):<15} "
            f"{r['total_length_gb']:<10} "
            f"{r['n50_mb']:<10} "
            f"{r['chromosomes_detected']:<10} "
            f"{expected:<10} "
            f"{match_str:<8}"
        )


def main() -> None:
    """Run batch processing example."""
    # In a real scenario, you would point to your assembly directory
    assembly_dir = Path("./assemblies")

    if not assembly_dir.exists():
        print("Example: Creating mock assemblies for demonstration")
        assembly_dir.mkdir(exist_ok=True)
        create_mock_assemblies(assembly_dir)

    results = batch_process(assembly_dir)

    if results:
        print_summary(results)

        # Save to JSON
        output_file = Path("batch_results.json")
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to {output_file}")


def create_mock_assemblies(assembly_dir: Path) -> None:
    """Create mock assembly files for demonstration."""
    for _name, info in SPECIES_INFO.items():
        fasta_path = assembly_dir / info["file"]
        n_chr = info["karyotype"]

        with open(fasta_path, "w") as f:
            # Create chromosome-sized scaffolds
            for i in range(1, n_chr + 1):
                f.write(f">chr{i}\n")
                size = 50_000_000 + (n_chr - i) * 5_000_000
                # Write in chunks to avoid huge strings
                remaining = size
                while remaining > 0:
                    chunk = min(80, remaining)
                    f.write("A" * chunk + "\n")
                    remaining -= chunk

            # Add some unplaced scaffolds
            for i in range(5):
                f.write(f">scaffold_unplaced_{i}\n")
                f.write("G" * 100_000 + "\n")


if __name__ == "__main__":
    main()
