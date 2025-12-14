"""
NCBI Assembly Report Parser for ChromDetect.

This module parses NCBI assembly report files to extract scaffold-chromosome
mappings, enabling more accurate scaffold classification.

NCBI assembly reports can be downloaded from:
  https://www.ncbi.nlm.nih.gov/assembly/

Example report header format:
  # Sequence-Name  Sequence-Role  Assigned-Molecule  ...  GenBank-Accn  ...
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class AssemblyReportEntry:
    """A single entry from an NCBI assembly report.

    Attributes:
        sequence_name: Name of the sequence (scaffold/contig)
        sequence_role: Role of sequence (assembled-molecule, unlocalized-scaffold, etc.)
        assigned_molecule: Chromosome assignment (1, 2, X, MT, na, etc.)
        assigned_molecule_type: Type of molecule (Chromosome, Mitochondrion, na)
        genbank_accession: GenBank accession number
        refseq_accession: RefSeq accession number
        length: Sequence length in bp (if available)
    """

    sequence_name: str
    sequence_role: str
    assigned_molecule: str
    assigned_molecule_type: str
    genbank_accession: str
    refseq_accession: str
    length: int | None = None


@dataclass
class AssemblyReport:
    """Parsed NCBI assembly report.

    Attributes:
        assembly_name: Name of the assembly
        organism: Organism name
        taxid: NCBI taxonomy ID
        entries: List of AssemblyReportEntry objects
        chromosome_map: Mapping of scaffold names to chromosome IDs
    """

    assembly_name: str | None
    organism: str | None
    taxid: str | None
    entries: list[AssemblyReportEntry]

    @property
    def chromosome_map(self) -> dict[str, str]:
        """Get mapping of scaffold names to chromosome IDs."""
        mapping = {}
        for entry in self.entries:
            if entry.assigned_molecule and entry.assigned_molecule.lower() != "na":
                # Map both sequence name and accessions
                mapping[entry.sequence_name] = entry.assigned_molecule
                if entry.genbank_accession and entry.genbank_accession.lower() != "na":
                    mapping[entry.genbank_accession] = entry.assigned_molecule
                if entry.refseq_accession and entry.refseq_accession.lower() != "na":
                    mapping[entry.refseq_accession] = entry.assigned_molecule
        return mapping

    @property
    def chromosome_scaffolds(self) -> set[str]:
        """Get set of scaffold names that are chromosome-level."""
        names = set()
        for entry in self.entries:
            if entry.sequence_role == "assembled-molecule":
                names.add(entry.sequence_name)
                if entry.genbank_accession and entry.genbank_accession.lower() != "na":
                    names.add(entry.genbank_accession)
                if entry.refseq_accession and entry.refseq_accession.lower() != "na":
                    names.add(entry.refseq_accession)
        return names

    @property
    def unlocalized_scaffolds(self) -> set[str]:
        """Get set of scaffold names that are unlocalized."""
        names = set()
        for entry in self.entries:
            if entry.sequence_role == "unlocalized-scaffold":
                names.add(entry.sequence_name)
                if entry.genbank_accession and entry.genbank_accession.lower() != "na":
                    names.add(entry.genbank_accession)
                if entry.refseq_accession and entry.refseq_accession.lower() != "na":
                    names.add(entry.refseq_accession)
        return names

    @property
    def unplaced_scaffolds(self) -> set[str]:
        """Get set of scaffold names that are unplaced."""
        names = set()
        for entry in self.entries:
            if entry.sequence_role == "unplaced-scaffold":
                names.add(entry.sequence_name)
                if entry.genbank_accession and entry.genbank_accession.lower() != "na":
                    names.add(entry.genbank_accession)
                if entry.refseq_accession and entry.refseq_accession.lower() != "na":
                    names.add(entry.refseq_accession)
        return names

    def get_expected_chromosome_count(self) -> int:
        """Get the expected number of chromosomes."""
        # Count unique assigned molecules that are chromosomes
        chromosomes = set()
        for entry in self.entries:
            if (
                entry.sequence_role == "assembled-molecule"
                and entry.assigned_molecule_type.lower() == "chromosome"
            ):
                chromosomes.add(entry.assigned_molecule)
        return len(chromosomes)


def parse_assembly_report(file_path: Path | str) -> AssemblyReport:
    """
    Parse an NCBI assembly report file.

    NCBI assembly reports are tab-separated files with a header that starts
    with '#'. The format includes columns like:
    - Sequence-Name
    - Sequence-Role
    - Assigned-Molecule
    - Assigned-Molecule-Location/Type
    - GenBank-Accn
    - Relationship
    - RefSeq-Accn
    - Assembly-Unit
    - Sequence-Length
    - UCSC-style-name

    Args:
        file_path: Path to the assembly report file

    Returns:
        AssemblyReport object with parsed data

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Assembly report file not found: {file_path}")

    assembly_name = None
    organism = None
    taxid = None
    entries = []
    column_indices: dict[str, int] = {}

    with open(file_path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n\r")

            if not line:
                continue

            # Parse comment/header lines
            if line.startswith("#"):
                line_content = line[1:].strip() if line.startswith("# ") else line[1:].strip()

                # Check if this is the header line (contains tab-separated column names)
                if "\t" in line_content and "sequence" in line_content.lower():
                    # Header line defines column positions
                    # Format: # Sequence-Name\tSequence-Role\t...
                    columns = line_content.split("\t")
                    for i, col in enumerate(columns):
                        col_name = col.strip().lower().replace("-", "_")
                        column_indices[col_name] = i
                # Check for metadata lines like "# Assembly name: GRCh38.p14"
                elif ":" in line_content:
                    key, _, value = line_content.partition(":")
                    key = key.strip().lower()
                    value = value.strip()
                    if key == "assembly name":
                        assembly_name = value
                    elif key == "organism name":
                        organism = value
                    elif key == "taxid":
                        taxid = value
                continue

            # Parse data line
            fields = line.split("\t")
            if len(fields) < 5:
                continue  # Skip malformed lines

            # Helper to get field by column name
            def get_field(
                name: str, field_list: list[str], default: str = ""
            ) -> str:
                idx = column_indices.get(name)
                if idx is not None and idx < len(field_list):
                    return field_list[idx].strip()
                return default

            # Parse length if available
            length_str = get_field("sequence_length", fields)
            length = None
            if length_str and length_str.lower() != "na":
                try:
                    length = int(length_str)
                except ValueError:
                    pass

            entry = AssemblyReportEntry(
                sequence_name=get_field("sequence_name", fields),
                sequence_role=get_field("sequence_role", fields),
                assigned_molecule=get_field("assigned_molecule", fields),
                assigned_molecule_type=get_field(
                    "assigned_molecule_location/type", fields
                )
                or get_field("assigned_molecule_type", fields),
                genbank_accession=get_field("genbank_accn", fields),
                refseq_accession=get_field("refseq_accn", fields),
                length=length,
            )

            if entry.sequence_name:
                entries.append(entry)

    if not entries:
        raise ValueError("No sequence entries found in assembly report")

    return AssemblyReport(
        assembly_name=assembly_name,
        organism=organism,
        taxid=taxid,
        entries=entries,
    )


def apply_assembly_report(
    scaffolds: list[tuple[str, int, str]],
    report: AssemblyReport,
) -> tuple[
    list[tuple[str, str, str | None]],  # (name, classification, chromosome_id)
    int,  # expected chromosome count
]:
    """
    Apply assembly report information to classify scaffolds.

    Args:
        scaffolds: List of (name, length, sequence) tuples from parse_fasta()
        report: Parsed AssemblyReport

    Returns:
        Tuple of (classifications, expected_chromosome_count)
        classifications is a list of (name, classification, chromosome_id) tuples
    """
    chr_scaffolds = report.chromosome_scaffolds
    unloc_scaffolds = report.unlocalized_scaffolds
    unplaced_scaffolds = report.unplaced_scaffolds
    chr_map = report.chromosome_map

    classifications = []
    for name, _length, _seq in scaffolds:
        if name in chr_scaffolds:
            chr_id = chr_map.get(name)
            classifications.append((name, "chromosome", chr_id))
        elif name in unloc_scaffolds:
            chr_id = chr_map.get(name)
            classifications.append((name, "unlocalized", chr_id))
        elif name in unplaced_scaffolds:
            classifications.append((name, "unplaced", None))
        else:
            # Not found in report
            classifications.append((name, "unknown", None))

    return classifications, report.get_expected_chromosome_count()
