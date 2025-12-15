"""
Core detection logic for ChromDetect.

This module contains the main classification algorithms for detecting
chromosome-level scaffolds in genome assemblies.
"""

from __future__ import annotations

import gzip
import typing
from dataclasses import asdict, dataclass
from pathlib import Path

from chromdetect.patterns import (
    COMPILED_CHROMOSOME_PATTERNS,
    COMPILED_FRAGMENT,
    COMPILED_UNLOCALIZED,
)


@dataclass
class ScaffoldInfo:
    """Information about a single scaffold.

    Attributes:
        name: Scaffold name from the FASTA header
        length: Scaffold length in base pairs
        classification: One of "chromosome", "unlocalized", "unplaced", "other"
        confidence: Confidence score from 0.0 to 1.0
        detection_method: Description of how classification was determined
        chromosome_id: Inferred chromosome ID if detected (e.g., "1", "X", "MT")
        gc_content: GC content as fraction (0.0-1.0), or None if not calculated
    """

    name: str
    length: int
    classification: str  # "chromosome", "unlocalized", "unplaced", "other"
    confidence: float  # 0.0 - 1.0
    detection_method: str
    chromosome_id: str | None = None
    gc_content: float | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


@dataclass
class AssemblyStats:
    """Summary statistics for the assembly.

    Attributes:
        total_scaffolds: Total number of scaffolds in assembly
        total_length: Total assembly length in base pairs
        n50: N50 scaffold length
        n90: N90 scaffold length
        chromosome_count: Number of scaffolds classified as chromosomes
        chromosome_length: Total length of chromosome-level scaffolds
        chromosome_n50: N50 of chromosome-level scaffolds only
        unlocalized_count: Number of unlocalized scaffolds
        unplaced_count: Number of unplaced scaffolds
        largest_scaffold: Length of largest scaffold
        gc_content: GC content as fraction (0.0-1.0) if calculated
    """

    total_scaffolds: int
    total_length: int
    n50: int
    n90: int
    chromosome_count: int
    chromosome_length: int
    chromosome_n50: int
    unlocalized_count: int
    unplaced_count: int
    largest_scaffold: int
    gc_content: float | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


def parse_fasta_from_handle(
    handle: typing.Iterable[str],
    keep_full_sequence: bool = False,
) -> list[tuple[str, int, str]]:
    """
    Parse FASTA from a file handle and return list of (name, length, sequence_sample).

    For efficiency, only stores first 10kb of each sequence (for GC calculation)
    unless keep_full_sequence is True.

    Args:
        handle: File-like object to read from (e.g., open file, sys.stdin)
        keep_full_sequence: If True, store entire sequence instead of just sample

    Returns:
        List of tuples: (scaffold_name, length, sequence or sequence_sample)

    Raises:
        ValueError: If the input appears to be empty or invalid FASTA format
    """
    scaffolds: list[tuple[str, int, str]] = []
    current_name: str | None = None
    current_length = 0
    current_seq_parts: list[str] = []
    sample_limit = 10000  # Only keep first 10kb for GC calculation (if not keeping full)
    line_count = 0
    found_header = False

    for raw_line in handle:
        line_count += 1
        line: str = str(raw_line).strip()

        # Skip empty lines
        if not line:
            continue

        if line.startswith(">"):
            found_header = True
            # Save previous scaffold
            if current_name is not None:
                seq = "".join(current_seq_parts)
                scaffolds.append((current_name, current_length, seq))

            # Start new scaffold - extract name (first word after >)
            header_parts = line[1:].split()
            if not header_parts:
                raise ValueError(
                    f"Line {line_count}: Empty FASTA header (just '>' with no name)"
                )
            current_name = header_parts[0]
            current_length = 0
            current_seq_parts = []
        else:
            # This is sequence data
            if not found_header:
                raise ValueError(
                    f"Line {line_count}: Sequence data before first FASTA header. "
                    "File may not be in FASTA format."
                )
            # Validate sequence characters (allow standard nucleotides + ambiguity codes)
            valid_chars = set("ACGTNUacgtnuRYSWKMBDHVryswkmbdhv.-")
            invalid_chars = set(line) - valid_chars
            if invalid_chars:
                # Only warn for non-whitespace invalid chars, don't fail
                # Some FASTA files have quality scores or other data
                pass
            current_length += len(line)
            if keep_full_sequence:
                current_seq_parts.append(line)
            elif sum(len(s) for s in current_seq_parts) < sample_limit:
                current_seq_parts.append(line)

    # Don't forget last scaffold
    if current_name is not None:
        seq = "".join(current_seq_parts)
        scaffolds.append((current_name, current_length, seq))

    if not scaffolds:
        if line_count == 0:
            raise ValueError("Input is empty")
        raise ValueError(
            "No scaffolds found in input. "
            "File may not be in FASTA format (expected lines starting with '>')."
        )

    return scaffolds


def parse_fasta(
    fasta_path: Path | str,
    keep_full_sequence: bool = False,
) -> list[tuple[str, int, str]]:
    """
    Parse FASTA file and return list of (name, length, sequence_sample).

    Handles gzipped files automatically based on .gz extension.
    For efficiency, only stores first 10kb of each sequence (for GC calculation)
    unless keep_full_sequence is True.

    Args:
        fasta_path: Path to FASTA file (can be gzipped)
        keep_full_sequence: If True, store entire sequence instead of just sample

    Returns:
        List of tuples: (scaffold_name, length, sequence or sequence_sample)

    Raises:
        FileNotFoundError: If the FASTA file doesn't exist
        ValueError: If the file appears to be empty or invalid
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    mode = "rt" if str(fasta_path).endswith(".gz") else "r"

    with opener(fasta_path, mode, encoding="utf-8") as f:
        scaffolds = parse_fasta_from_handle(f, keep_full_sequence)  # type: ignore[arg-type]

    return scaffolds


def calculate_gc(sequence: str) -> float | None:
    """Calculate GC content of a sequence.

    Args:
        sequence: DNA sequence string

    Returns:
        GC content as fraction (0.0-1.0), or None if sequence is empty
    """
    if not sequence:
        return None

    sequence = sequence.upper()
    gc = sequence.count("G") + sequence.count("C")
    total = gc + sequence.count("A") + sequence.count("T")

    return gc / total if total > 0 else None


def calculate_n50(lengths: list[int]) -> int:
    """Calculate N50 from list of scaffold lengths.

    N50 is the length such that 50% of the total assembly length
    is contained in scaffolds of this length or longer.

    Args:
        lengths: List of scaffold lengths

    Returns:
        N50 value, or 0 if list is empty
    """
    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    running_sum = 0

    for length in sorted_lengths:
        running_sum += length
        if running_sum >= total / 2:
            return length

    return sorted_lengths[-1]


def calculate_n90(lengths: list[int]) -> int:
    """Calculate N90 from list of scaffold lengths.

    Args:
        lengths: List of scaffold lengths

    Returns:
        N90 value, or 0 if list is empty
    """
    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    running_sum = 0

    for length in sorted_lengths:
        running_sum += length
        if running_sum >= total * 0.9:
            return length

    return sorted_lengths[-1]


def detect_by_name(
    name: str,
    custom_patterns: tuple | None = None,
) -> tuple[str, float, str, str | None]:
    """
    Detect scaffold type by name pattern matching.

    Uses pre-compiled regex patterns to identify chromosome-level scaffolds
    based on common naming conventions.

    Args:
        name: Scaffold name from FASTA header
        custom_patterns: Optional tuple of (chr_patterns, unloc_patterns, frag_patterns)
                        from merge_patterns(). If None, uses built-in patterns.

    Returns:
        Tuple of (classification, confidence, method, chromosome_id)
        - classification: "chromosome", "unlocalized", "unplaced", or "other"
        - confidence: 0.0-1.0 confidence score
        - method: Description of detection method used
        - chromosome_id: Extracted chromosome ID if available (e.g., "1", "X")
    """
    # Use custom patterns if provided, otherwise use built-in
    if custom_patterns:
        chr_patterns, unloc_patterns, frag_patterns = custom_patterns
    else:
        chr_patterns = COMPILED_CHROMOSOME_PATTERNS
        unloc_patterns = COMPILED_UNLOCALIZED
        frag_patterns = COMPILED_FRAGMENT

    # Check for unlocalized patterns first (these override chromosome patterns)
    for pattern in unloc_patterns:
        if pattern.search(name):
            return ("unlocalized", 0.8, "name_unlocalized", None)

    # Check for fragment patterns
    for pattern in frag_patterns:
        if pattern.search(name):
            return ("unplaced", 0.6, "name_fragment", None)

    # Check for chromosome patterns
    for pattern, method in chr_patterns:
        match = pattern.match(name)
        if match:
            chr_id = match.group(1) if match.lastindex else None
            return ("chromosome", 0.9, f"name_{method}", chr_id)

    # Default: unknown
    return ("other", 0.3, "name_none", None)


def detect_by_size(
    length: int,
    n50: int,
    largest: int,
    min_chromosome_size: int = 10_000_000,
) -> tuple[str, float, str]:
    """
    Detect scaffold type by size heuristics.

    Large scaffolds are typically chromosomes. This function uses size relative
    to N50 and assembly-specific thresholds.

    Args:
        length: Scaffold length in base pairs
        n50: Assembly N50 value
        largest: Length of largest scaffold
        min_chromosome_size: Minimum size to consider chromosome-level (default 10Mb)

    Returns:
        Tuple of (classification, confidence, method)
    """
    # Very large scaffolds are almost certainly chromosomes
    if length >= min_chromosome_size:
        # Scale confidence by how close to largest
        confidence = 0.7 + 0.3 * (length / largest) if largest > 0 else 0.7
        return ("chromosome", min(confidence, 0.95), "size_large")

    # Scaffolds near N50 are likely chromosomes
    if length >= n50 * 0.5:
        confidence = 0.5 + 0.3 * (length / n50) if n50 > 0 else 0.5
        return ("chromosome", min(confidence, 0.8), "size_n50")

    # Small scaffolds are likely unplaced
    if length < 1_000_000:
        return ("unplaced", 0.6, "size_small")

    # Medium scaffolds are ambiguous
    return ("other", 0.4, "size_medium")


def classify_scaffolds(
    scaffolds: list[tuple[str, int, str]],
    min_chromosome_size: int = 10_000_000,
    expected_chromosomes: int | None = None,
    custom_patterns: tuple | None = None,
    assembly_report: typing.Any | None = None,
) -> tuple[list[ScaffoldInfo], AssemblyStats]:
    """
    Classify all scaffolds using multiple detection strategies.

    This is the main entry point for scaffold classification. It combines
    name-based and size-based detection with optional karyotype adjustment.

    Args:
        scaffolds: List of (name, length, sequence_sample) tuples from parse_fasta()
        min_chromosome_size: Minimum size in bp to consider chromosome-level
        expected_chromosomes: Known chromosome count (karyotype) for adjustment
        custom_patterns: Optional tuple of (chr_patterns, unloc_patterns, frag_patterns)
                        from merge_patterns(). If None, uses built-in patterns.
        assembly_report: Optional AssemblyReport object for NCBI-based classification.
                        If provided, uses the report for authoritative classification.

    Returns:
        Tuple of (list of ScaffoldInfo, AssemblyStats)

    Raises:
        ValueError: If no scaffolds provided

    Example:
        >>> scaffolds = parse_fasta("assembly.fasta")
        >>> results, stats = classify_scaffolds(scaffolds)
        >>> print(f"Found {stats.chromosome_count} chromosomes")
    """
    if not scaffolds:
        raise ValueError("No scaffolds found in assembly")

    lengths = [s[1] for s in scaffolds]
    n50 = calculate_n50(lengths)
    n90 = calculate_n90(lengths)
    largest = max(lengths)
    total_length = sum(lengths)

    # Build scaffold sequence lookup for GC calculation
    scaffold_seqs = {name: seq for name, _length, seq in scaffolds}

    # Pre-compute assembly report classifications if provided
    report_classifications: dict[str, tuple[str, str | None]] = {}
    if assembly_report is not None:
        # assembly_report is an AssemblyReport object
        from chromdetect.assembly_report import apply_assembly_report

        report_class_list, report_expected = apply_assembly_report(
            scaffolds, assembly_report
        )
        for name, classification, chr_id in report_class_list:
            if classification != "unknown":
                report_classifications[name] = (classification, chr_id)
        # Use report's expected count if not already specified
        if expected_chromosomes is None and report_expected > 0:
            expected_chromosomes = report_expected

    results = []

    for name, length, _seq_sample in scaffolds:
        # Check if we have an authoritative classification from assembly report
        if name in report_classifications:
            report_class, report_chr_id = report_classifications[name]
            # Assembly report is authoritative - use it directly
            final_class = report_class
            final_conf = 0.99  # Very high confidence from NCBI
            final_method = "ncbi_report"
            chr_id = report_chr_id

            # Calculate GC content for this scaffold
            scaffold_gc = None
            if name in scaffold_seqs:
                seq = scaffold_seqs[name]
                if seq:
                    scaffold_gc = calculate_gc(seq)
                    if scaffold_gc is not None:
                        scaffold_gc = round(scaffold_gc, 4)

            results.append(
                ScaffoldInfo(
                    name=name,
                    length=length,
                    classification=final_class,
                    confidence=round(final_conf, 3),
                    detection_method=final_method,
                    chromosome_id=chr_id,
                    gc_content=scaffold_gc,
                )
            )
            continue

        # Get classifications from each method
        name_class, name_conf, name_method, chr_id = detect_by_name(
            name, custom_patterns
        )
        size_class, size_conf, size_method = detect_by_size(
            length, n50, largest, min_chromosome_size
        )

        # Combine classifications with priority rules
        if name_conf >= 0.8:
            # Strong name-based classification takes priority
            final_class = name_class
            final_conf = name_conf
            final_method = name_method
        elif size_conf >= 0.7 and size_class == "chromosome":
            # Strong size-based chromosome detection
            if name_class == "chromosome":
                # Both methods agree - boost confidence
                final_class = "chromosome"
                final_conf = min(0.95, (name_conf + size_conf) / 2 + 0.1)
                final_method = f"{name_method}+{size_method}"
            else:
                # Size says chromosome, name doesn't - use size with penalty
                final_class = "chromosome"
                final_conf = size_conf * 0.9
                final_method = size_method
        elif name_class == "chromosome":
            # Weak name-based chromosome
            final_class = "chromosome"
            final_conf = name_conf
            final_method = name_method
        elif size_class == "unplaced":
            final_class = "unplaced"
            final_conf = size_conf
            final_method = size_method
        else:
            # Default to name-based with low confidence
            final_class = name_class if name_class != "other" else "unplaced"
            final_conf = max(name_conf, size_conf) * 0.8
            final_method = name_method

        # Calculate GC content for this scaffold
        scaffold_gc = None
        if name in scaffold_seqs:
            seq = scaffold_seqs[name]
            if seq:
                scaffold_gc = calculate_gc(seq)
                if scaffold_gc is not None:
                    scaffold_gc = round(scaffold_gc, 4)

        results.append(
            ScaffoldInfo(
                name=name,
                length=length,
                classification=final_class,
                confidence=round(final_conf, 3),
                detection_method=final_method,
                chromosome_id=chr_id,
                gc_content=scaffold_gc,
            )
        )

    # Adjust for expected karyotype if provided
    if expected_chromosomes is not None:
        results = _adjust_for_karyotype(results, expected_chromosomes)

    # Calculate statistics
    chromosomes = [r for r in results if r.classification == "chromosome"]
    unlocalized = [r for r in results if r.classification == "unlocalized"]
    chr_lengths = [r.length for r in chromosomes]

    # Calculate overall GC from samples
    all_seqs = "".join(s[2] for s in scaffolds[:100])
    gc_content = calculate_gc(all_seqs)

    stats = AssemblyStats(
        total_scaffolds=len(results),
        total_length=total_length,
        n50=n50,
        n90=n90,
        chromosome_count=len(chromosomes),
        chromosome_length=sum(chr_lengths),
        chromosome_n50=calculate_n50(chr_lengths) if chr_lengths else 0,
        unlocalized_count=len(unlocalized),
        unplaced_count=len([r for r in results if r.classification == "unplaced"]),
        largest_scaffold=largest,
        gc_content=round(gc_content, 4) if gc_content else None,
    )

    return results, stats


def _adjust_for_karyotype(
    results: list[ScaffoldInfo],
    expected: int,
) -> list[ScaffoldInfo]:
    """
    Adjust classifications based on expected chromosome count.

    If we have more chromosome candidates than expected, demote lowest confidence.
    If we have fewer, promote largest unplaced scaffolds.

    Args:
        results: List of ScaffoldInfo from initial classification
        expected: Expected chromosome count from karyotype

    Returns:
        Adjusted list of ScaffoldInfo
    """
    chromosomes = [r for r in results if r.classification == "chromosome"]
    current_count = len(chromosomes)

    if current_count == expected:
        return results

    # Sort all results for potential adjustment
    results_sorted = sorted(results, key=lambda r: (-r.length, -r.confidence))

    if current_count > expected:
        # Too many chromosomes - demote lowest confidence ones
        chr_by_conf = sorted(chromosomes, key=lambda r: r.confidence)
        to_demote = current_count - expected

        demote_names = {r.name for r in chr_by_conf[:to_demote]}

        for r in results:
            if r.name in demote_names:
                r.classification = "unplaced"
                r.detection_method += "_demoted_karyotype"
                r.confidence *= 0.5

    else:
        # Too few chromosomes - promote largest unplaced
        unplaced = [r for r in results_sorted if r.classification == "unplaced"]
        to_promote = expected - current_count

        promote_names = {r.name for r in unplaced[:to_promote]}

        for r in results:
            if r.name in promote_names:
                r.classification = "chromosome"
                r.detection_method += "_promoted_karyotype"
                r.confidence = min(0.6, r.confidence + 0.2)

    return results


def write_fasta(
    sequences: list[tuple[str, str]],
    output_path: Path | str | None = None,
    line_width: int = 80,
) -> str:
    """
    Write sequences in FASTA format.

    Args:
        sequences: List of (name, sequence) tuples
        output_path: Path to write to, or None to return as string
        line_width: Characters per line for sequence wrapping (default 80)

    Returns:
        FASTA-formatted string if output_path is None, otherwise empty string
    """
    lines = []
    for name, seq in sequences:
        lines.append(f">{name}")
        # Wrap sequence at specified line width
        for i in range(0, len(seq), line_width):
            lines.append(seq[i : i + line_width])

    fasta_content = "\n".join(lines)
    if fasta_content:
        fasta_content += "\n"

    if output_path:
        output_path = Path(output_path)
        with open(output_path, "w") as f:
            f.write(fasta_content)
        return ""

    return fasta_content


def format_bed(
    results: list[ScaffoldInfo],
    include_header: bool = False,
) -> str:
    """
    Format scaffold information as BED format.

    BED format is 0-based, half-open coordinates. Each scaffold becomes
    a region spanning its full length.

    Args:
        results: List of ScaffoldInfo from classification
        include_header: If True, include a header line (not standard BED)

    Returns:
        BED-formatted string
    """
    lines = []
    if include_header:
        lines.append("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand")

    for r in results:
        # BED format: chrom, chromStart (0-based), chromEnd, name, score, strand
        # Use confidence * 1000 for score (BED scores are 0-1000)
        score = int(r.confidence * 1000)
        # Use classification as the name field, scaffold name as chrom
        lines.append(f"{r.name}\t0\t{r.length}\t{r.classification}\t{score}\t.")

    return "\n".join(lines)


def format_gff(
    results: list[ScaffoldInfo],
    source: str = "chromdetect",
) -> str:
    """
    Format scaffold information as GFF3 format.

    Each scaffold becomes a region feature spanning its full length.

    Args:
        results: List of ScaffoldInfo from classification
        source: Source field for GFF (default "chromdetect")

    Returns:
        GFF3-formatted string
    """
    lines = ["##gff-version 3"]

    for r in results:
        # GFF3 format: seqid, source, type, start (1-based), end, score, strand, phase, attributes
        feature_type = "chromosome" if r.classification == "chromosome" else "scaffold"
        # GFF scores are . or float
        score = f"{r.confidence:.3f}"
        # Attributes: key=value pairs separated by semicolons
        attrs = [
            f"ID={r.name}",
            f"Name={r.name}",
            f"classification={r.classification}",
            f"detection_method={r.detection_method}",
        ]
        if r.chromosome_id:
            attrs.append(f"chromosome_id={r.chromosome_id}")

        attr_str = ";".join(attrs)
        lines.append(
            f"{r.name}\t{source}\t{feature_type}\t1\t{r.length}\t{score}\t.\t.\t{attr_str}"
        )

    return "\n".join(lines)


def classify_fasta(
    fasta_path: str | Path,
    min_chromosome_size: int = 10_000_000,
    expected_chromosomes: int | None = None,
) -> tuple[list[ScaffoldInfo], AssemblyStats]:
    """
    Convenience function to classify scaffolds directly from a FASTA file.

    This is a high-level API that combines parse_fasta() and classify_scaffolds()
    into a single call for simpler usage.

    Args:
        fasta_path: Path to FASTA file (can be gzipped)
        min_chromosome_size: Minimum size in bp to consider chromosome-level
        expected_chromosomes: Known chromosome count (karyotype) for adjustment

    Returns:
        Tuple of (list of ScaffoldInfo, AssemblyStats)

    Example:
        >>> results, stats = classify_fasta("assembly.fasta")
        >>> print(f"Found {stats.chromosome_count} chromosomes")
    """
    scaffolds = parse_fasta(fasta_path)

    return classify_scaffolds(
        scaffolds,
        min_chromosome_size=min_chromosome_size,
        expected_chromosomes=expected_chromosomes,
    )
