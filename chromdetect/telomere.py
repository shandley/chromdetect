"""
Telomere detection for ChromDetect.

This module provides functions to detect telomeric repeats at scaffold ends,
which is a strong indicator of chromosome-level assembly.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Common telomere repeat motifs across different organisms
# Format: (forward_motif, reverse_complement, description)
TELOMERE_MOTIFS = [
    # Vertebrates, fungi, and most eukaryotes
    ("TTAGGG", "CCCTAA", "vertebrate"),
    # Plants (Arabidopsis-type)
    ("TTTAGGG", "CCCTAAA", "plant_arabidopsis"),
    # Some insects (Bombyx mori)
    ("TTAGG", "CCTAA", "insect_bombyx"),
    # Nematodes (C. elegans)
    ("TTAGGC", "GCCTAA", "nematode"),
    # Some ciliates
    ("TTGGGG", "CCCCAA", "ciliate_tetrahymena"),
    ("TTTTGGGG", "CCCCAAAA", "ciliate_oxytricha"),
    # Plasmodium
    ("TTTTAGGG", "CCCTAAAA", "plasmodium"),
    # Some green algae
    ("TTTTAGGG", "CCCTAAAA", "green_algae"),
]

# Minimum number of consecutive repeats to consider as telomeric
MIN_REPEATS = 3

# How many bp from scaffold ends to search for telomeres
DEFAULT_SEARCH_WINDOW = 10000


@dataclass
class TelomereResult:
    """Result of telomere detection for a scaffold.

    Attributes:
        has_5prime: True if telomere detected at 5' end
        has_3prime: True if telomere detected at 3' end
        motif_5prime: Detected motif at 5' end (if any)
        motif_3prime: Detected motif at 3' end (if any)
        repeats_5prime: Number of repeats at 5' end
        repeats_3prime: Number of repeats at 3' end
        organism_type: Inferred organism type based on motif
    """

    has_5prime: bool = False
    has_3prime: bool = False
    motif_5prime: str | None = None
    motif_3prime: str | None = None
    repeats_5prime: int = 0
    repeats_3prime: int = 0
    organism_type: str | None = None

    @property
    def has_telomere(self) -> bool:
        """Return True if telomere detected at either end."""
        return self.has_5prime or self.has_3prime

    @property
    def is_complete(self) -> bool:
        """Return True if telomeres detected at both ends (T2T)."""
        return self.has_5prime and self.has_3prime

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "has_5prime": self.has_5prime,
            "has_3prime": self.has_3prime,
            "motif_5prime": self.motif_5prime,
            "motif_3prime": self.motif_3prime,
            "repeats_5prime": self.repeats_5prime,
            "repeats_3prime": self.repeats_3prime,
            "organism_type": self.organism_type,
            "has_telomere": self.has_telomere,
            "is_complete": self.is_complete,
        }


def _count_repeats(sequence: str, motif: str) -> int:
    """Count consecutive repeats of a motif in a sequence.

    Args:
        sequence: DNA sequence to search
        motif: Telomere motif to count

    Returns:
        Maximum number of consecutive repeats found
    """
    if not sequence or not motif:
        return 0

    sequence = sequence.upper()
    motif = motif.upper()

    # Build regex for consecutive repeats
    # Match the motif repeated 1 or more times
    pattern = f"({motif})+"
    matches = re.finditer(pattern, sequence, re.IGNORECASE)

    max_repeats = 0
    for match in matches:
        match_len = len(match.group())
        repeats = match_len // len(motif)
        max_repeats = max(max_repeats, repeats)

    return max_repeats


def detect_telomere(
    sequence: str,
    search_window: int = DEFAULT_SEARCH_WINDOW,
    min_repeats: int = MIN_REPEATS,
    motifs: list[tuple[str, str, str]] | None = None,
) -> TelomereResult:
    """Detect telomeric repeats at scaffold ends.

    Searches both ends of the sequence for known telomere motifs.
    At the 5' end, looks for reverse complement (CCCTAA for vertebrates).
    At the 3' end, looks for forward motif (TTAGGG for vertebrates).

    Args:
        sequence: Full scaffold sequence
        search_window: Number of bp from each end to search (default 10kb)
        min_repeats: Minimum consecutive repeats to count as telomere
        motifs: Custom motifs list, or None to use defaults

    Returns:
        TelomereResult with detection information
    """
    if not sequence:
        return TelomereResult()

    sequence = sequence.upper()
    seq_len = len(sequence)

    # Get sequence ends
    window = min(search_window, seq_len // 2)
    five_prime = sequence[:window]
    three_prime = sequence[-window:] if seq_len > window else sequence

    if motifs is None:
        motifs = TELOMERE_MOTIFS

    result = TelomereResult()
    best_5prime_repeats = 0
    best_3prime_repeats = 0

    for forward, reverse, org_type in motifs:
        # 5' end typically has reverse complement (CCCTAA...)
        repeats_5prime = _count_repeats(five_prime, reverse)
        # 3' end typically has forward motif (...TTAGGG)
        repeats_3prime = _count_repeats(three_prime, forward)

        # Update if this motif has more repeats
        if repeats_5prime >= min_repeats and repeats_5prime > best_5prime_repeats:
            best_5prime_repeats = repeats_5prime
            result.has_5prime = True
            result.motif_5prime = reverse
            result.repeats_5prime = repeats_5prime
            result.organism_type = org_type

        if repeats_3prime >= min_repeats and repeats_3prime > best_3prime_repeats:
            best_3prime_repeats = repeats_3prime
            result.has_3prime = True
            result.motif_3prime = forward
            result.repeats_3prime = repeats_3prime
            # Only set organism type if not already set, or if this has more total
            if result.organism_type is None:
                result.organism_type = org_type

    return result


def detect_telomeres_batch(
    scaffolds: list[tuple[str, int, str]],
    search_window: int = DEFAULT_SEARCH_WINDOW,
    min_repeats: int = MIN_REPEATS,
) -> dict[str, TelomereResult]:
    """Detect telomeres for multiple scaffolds.

    Args:
        scaffolds: List of (name, length, sequence) tuples
        search_window: Number of bp from each end to search
        min_repeats: Minimum consecutive repeats to count as telomere

    Returns:
        Dictionary mapping scaffold name to TelomereResult
    """
    results = {}
    for name, _length, sequence in scaffolds:
        results[name] = detect_telomere(sequence, search_window, min_repeats)
    return results


def get_telomere_summary(telomere_results: dict[str, TelomereResult]) -> dict:
    """Generate summary statistics for telomere detection.

    Args:
        telomere_results: Dictionary from detect_telomeres_batch

    Returns:
        Summary dictionary with counts and percentages
    """
    total = len(telomere_results)
    if total == 0:
        return {
            "total_scaffolds": 0,
            "with_telomere": 0,
            "with_5prime": 0,
            "with_3prime": 0,
            "complete_t2t": 0,
            "percent_with_telomere": 0.0,
            "percent_complete": 0.0,
        }

    with_telomere = sum(1 for r in telomere_results.values() if r.has_telomere)
    with_5prime = sum(1 for r in telomere_results.values() if r.has_5prime)
    with_3prime = sum(1 for r in telomere_results.values() if r.has_3prime)
    complete = sum(1 for r in telomere_results.values() if r.is_complete)

    return {
        "total_scaffolds": total,
        "with_telomere": with_telomere,
        "with_5prime": with_5prime,
        "with_3prime": with_3prime,
        "complete_t2t": complete,
        "percent_with_telomere": round(100 * with_telomere / total, 1),
        "percent_complete": round(100 * complete / total, 1),
    }
