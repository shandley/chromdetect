"""
Tests for NCBI assembly report parsing functionality.
"""


import pytest

from chromdetect.assembly_report import (
    AssemblyReportEntry,
    apply_assembly_report,
    parse_assembly_report,
)


@pytest.fixture
def sample_report_content():
    """Sample NCBI assembly report content."""
    return """\
# Assembly name:  Test_Assembly_v1
# Organism name:  Test organism
# Taxid:          12345
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
chr1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	248956422	chr1
chr2	assembled-molecule	2	Chromosome	CM000002.1	=	NC_000002.1	Primary Assembly	242193529	chr2
chrX	assembled-molecule	X	Chromosome	CM000023.1	=	NC_000023.1	Primary Assembly	156040895	chrX
chrM	assembled-molecule	MT	Mitochondrion	CM000024.1	=	NC_012920.1	non-nuclear	16569	chrM
chr1_unlocalized_001	unlocalized-scaffold	1	Chromosome	CM000001.1	<>	na	Primary Assembly	50000	na
scaffold_001	unplaced-scaffold	na	na	CM000025.1	=	na	Primary Assembly	25000	na
"""


@pytest.fixture
def sample_report_file(tmp_path, sample_report_content):
    """Create a sample report file."""
    report_file = tmp_path / "assembly_report.txt"
    report_file.write_text(sample_report_content)
    return report_file


class TestAssemblyReportEntry:
    """Tests for AssemblyReportEntry dataclass."""

    def test_entry_creation(self):
        """Test creating an entry."""
        entry = AssemblyReportEntry(
            sequence_name="chr1",
            sequence_role="assembled-molecule",
            assigned_molecule="1",
            assigned_molecule_type="Chromosome",
            genbank_accession="CM000001.1",
            refseq_accession="NC_000001.1",
            length=248956422,
        )
        assert entry.sequence_name == "chr1"
        assert entry.assigned_molecule == "1"
        assert entry.length == 248956422


class TestParseAssemblyReport:
    """Tests for parsing NCBI assembly reports."""

    def test_parse_basic_report(self, sample_report_file):
        """Test parsing a basic assembly report."""
        report = parse_assembly_report(sample_report_file)

        assert report.assembly_name == "Test_Assembly_v1"
        assert report.organism == "Test organism"
        assert report.taxid == "12345"
        assert len(report.entries) == 6

    def test_parse_chromosome_entries(self, sample_report_file):
        """Test that chromosome entries are correctly parsed."""
        report = parse_assembly_report(sample_report_file)

        # Find chr1 entry
        chr1 = next(e for e in report.entries if e.sequence_name == "chr1")
        assert chr1.sequence_role == "assembled-molecule"
        assert chr1.assigned_molecule == "1"
        assert chr1.assigned_molecule_type == "Chromosome"
        assert chr1.genbank_accession == "CM000001.1"
        assert chr1.refseq_accession == "NC_000001.1"

    def test_parse_unlocalized_entries(self, sample_report_file):
        """Test that unlocalized entries are correctly parsed."""
        report = parse_assembly_report(sample_report_file)

        unloc = next(e for e in report.entries if "unlocalized" in e.sequence_name)
        assert unloc.sequence_role == "unlocalized-scaffold"
        assert unloc.assigned_molecule == "1"

    def test_parse_unplaced_entries(self, sample_report_file):
        """Test that unplaced entries are correctly parsed."""
        report = parse_assembly_report(sample_report_file)

        unplaced = next(e for e in report.entries if e.sequence_name == "scaffold_001")
        assert unplaced.sequence_role == "unplaced-scaffold"
        assert unplaced.assigned_molecule == "na"

    def test_file_not_found(self, tmp_path):
        """Test error for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_assembly_report(tmp_path / "nonexistent.txt")

    def test_empty_file(self, tmp_path):
        """Test error for empty file."""
        empty_file = tmp_path / "empty.txt"
        empty_file.write_text("")
        with pytest.raises(ValueError, match="No sequence entries"):
            parse_assembly_report(empty_file)


class TestAssemblyReportProperties:
    """Tests for AssemblyReport computed properties."""

    def test_chromosome_map(self, sample_report_file):
        """Test chromosome mapping."""
        report = parse_assembly_report(sample_report_file)
        chr_map = report.chromosome_map

        # Should map sequence names
        assert chr_map["chr1"] == "1"
        assert chr_map["chr2"] == "2"
        assert chr_map["chrX"] == "X"
        assert chr_map["chrM"] == "MT"

        # Should map accessions too
        assert chr_map["NC_000001.1"] == "1"
        assert chr_map["CM000001.1"] == "1"

        # Unlocalized should also be mapped
        assert chr_map["chr1_unlocalized_001"] == "1"

    def test_chromosome_scaffolds(self, sample_report_file):
        """Test chromosome scaffold set."""
        report = parse_assembly_report(sample_report_file)
        chr_scaffolds = report.chromosome_scaffolds

        assert "chr1" in chr_scaffolds
        assert "chr2" in chr_scaffolds
        assert "chrX" in chr_scaffolds
        assert "chrM" in chr_scaffolds
        assert "NC_000001.1" in chr_scaffolds

        # Unlocalized should not be in chromosome scaffolds
        assert "chr1_unlocalized_001" not in chr_scaffolds

    def test_unlocalized_scaffolds(self, sample_report_file):
        """Test unlocalized scaffold set."""
        report = parse_assembly_report(sample_report_file)
        unloc_scaffolds = report.unlocalized_scaffolds

        assert "chr1_unlocalized_001" in unloc_scaffolds
        assert "chr1" not in unloc_scaffolds

    def test_unplaced_scaffolds(self, sample_report_file):
        """Test unplaced scaffold set."""
        report = parse_assembly_report(sample_report_file)
        unplaced = report.unplaced_scaffolds

        assert "scaffold_001" in unplaced
        assert "chr1" not in unplaced

    def test_expected_chromosome_count(self, sample_report_file):
        """Test expected chromosome count."""
        report = parse_assembly_report(sample_report_file)

        # Should count unique chromosome assignments (1, 2, X)
        # MT is Mitochondrion, not Chromosome
        count = report.get_expected_chromosome_count()
        assert count == 3  # chr1, chr2, chrX


class TestApplyAssemblyReport:
    """Tests for applying assembly report to scaffolds."""

    def test_apply_to_scaffolds(self, sample_report_file):
        """Test applying report to scaffold list."""
        report = parse_assembly_report(sample_report_file)

        scaffolds = [
            ("chr1", 248956422, "ACGT" * 100),
            ("chr2", 242193529, "ACGT" * 100),
            ("scaffold_001", 25000, "ACGT" * 100),
            ("unknown_scaffold", 10000, "ACGT" * 100),
        ]

        classifications, expected = apply_assembly_report(scaffolds, report)

        assert len(classifications) == 4
        assert classifications[0] == ("chr1", "chromosome", "1")
        assert classifications[1] == ("chr2", "chromosome", "2")
        assert classifications[2] == ("scaffold_001", "unplaced", None)
        assert classifications[3] == ("unknown_scaffold", "unknown", None)

    def test_apply_with_accessions(self, sample_report_file):
        """Test applying report when scaffolds use accession names."""
        report = parse_assembly_report(sample_report_file)

        scaffolds = [
            ("NC_000001.1", 248956422, "ACGT" * 100),
            ("CM000002.1", 242193529, "ACGT" * 100),
        ]

        classifications, expected = apply_assembly_report(scaffolds, report)

        # Should still classify correctly using accession mapping
        assert classifications[0] == ("NC_000001.1", "chromosome", "1")
        assert classifications[1] == ("CM000002.1", "chromosome", "2")

    def test_expected_count(self, sample_report_file):
        """Test expected chromosome count from report."""
        report = parse_assembly_report(sample_report_file)

        scaffolds = [("chr1", 100000, "ACGT")]
        _, expected = apply_assembly_report(scaffolds, report)

        assert expected == 3  # chr1, chr2, chrX (not MT)
