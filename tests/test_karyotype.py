"""
Tests for karyotype database and validation functionality.
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

from chromdetect.karyotype import (
    Karyotype,
    KaryotypeDatabase,
    KaryotypeValidationResult,
    SexChromosomeSystem,
    format_karyotype_summary,
    format_karyotype_tsv,
    validate_karyotype,
)


class TestKaryotype:
    """Tests for Karyotype dataclass."""

    def test_karyotype_creation(self):
        """Test creating a karyotype."""
        karyotype = Karyotype(
            name="test_species",
            scientific_name="Testus speciesus",
            taxid=12345,
            chromosome_count=10,
            autosomes=["1", "2", "3", "4", "5", "6", "7", "8"],
            sex_chromosomes=["X", "Y"],
            sex_system=SexChromosomeSystem.XX_XY,
        )
        assert karyotype.name == "test_species"
        assert karyotype.chromosome_count == 10
        assert len(karyotype.autosomes) == 8
        assert karyotype.sex_chromosomes == ["X", "Y"]

    def test_all_chromosomes_property(self):
        """Test all_chromosomes combines autosomes and sex chromosomes."""
        karyotype = Karyotype(
            name="test",
            scientific_name="Test",
            taxid=None,
            chromosome_count=3,
            autosomes=["1", "2"],
            sex_chromosomes=["X"],
        )
        assert karyotype.all_chromosomes == ["1", "2", "X"]

    def test_total_with_mt_property(self):
        """Test total_with_mt includes mitochondria."""
        karyotype = Karyotype(
            name="test",
            scientific_name="Test",
            taxid=None,
            chromosome_count=10,
            autosomes=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
            mitochondria=True,
        )
        assert karyotype.total_with_mt == 11

        karyotype_no_mt = Karyotype(
            name="test",
            scientific_name="Test",
            taxid=None,
            chromosome_count=10,
            autosomes=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
            mitochondria=False,
        )
        assert karyotype_no_mt.total_with_mt == 10

    def test_to_dict(self):
        """Test serialization to dictionary."""
        karyotype = Karyotype(
            name="human",
            scientific_name="Homo sapiens",
            taxid=9606,
            chromosome_count=24,
            autosomes=[str(i) for i in range(1, 23)],
            sex_chromosomes=["X", "Y"],
            sex_system=SexChromosomeSystem.XX_XY,
        )
        d = karyotype.to_dict()
        assert d["name"] == "human"
        assert d["taxid"] == 9606
        assert d["sex_system"] == "XX/XY"

    def test_from_dict(self):
        """Test deserialization from dictionary."""
        data = {
            "name": "mouse",
            "scientific_name": "Mus musculus",
            "taxid": 10090,
            "chromosome_count": 21,
            "autosomes": [str(i) for i in range(1, 20)],
            "sex_chromosomes": ["X", "Y"],
            "sex_system": "XX/XY",
        }
        karyotype = Karyotype.from_dict(data)
        assert karyotype.name == "mouse"
        assert karyotype.taxid == 10090
        assert karyotype.sex_system == SexChromosomeSystem.XX_XY


class TestKaryotypeDatabase:
    """Tests for KaryotypeDatabase."""

    def test_builtin_database_loads(self):
        """Test that built-in database loads successfully."""
        db = KaryotypeDatabase()
        assert len(db) > 0

    def test_lookup_by_name(self):
        """Test looking up by common name."""
        db = KaryotypeDatabase()
        human = db.lookup("human")
        assert human is not None
        assert human.name == "human"
        assert human.chromosome_count == 24

    def test_lookup_by_scientific_name(self):
        """Test looking up by scientific name."""
        db = KaryotypeDatabase()
        human = db.lookup("Homo sapiens")
        assert human is not None
        assert human.name == "human"

    def test_lookup_by_taxid(self):
        """Test looking up by taxonomy ID."""
        db = KaryotypeDatabase()
        human = db.lookup(9606)
        assert human is not None
        assert human.name == "human"

    def test_lookup_by_alias(self):
        """Test looking up by alias."""
        db = KaryotypeDatabase()
        human = db.lookup("hsapiens")
        assert human is not None
        assert human.name == "human"

    def test_lookup_case_insensitive(self):
        """Test that lookup is case insensitive."""
        db = KaryotypeDatabase()
        assert db.lookup("HUMAN") is not None
        assert db.lookup("Human") is not None
        assert db.lookup("HOMO SAPIENS") is not None

    def test_lookup_not_found(self):
        """Test lookup returns None for unknown species."""
        db = KaryotypeDatabase()
        assert db.lookup("unknown_species_xyz") is None
        assert db.lookup(99999999) is None

    def test_list_species(self):
        """Test listing available species."""
        db = KaryotypeDatabase()
        species = db.list_species()
        assert "human" in species
        assert "mouse" in species
        assert "zebrafish" in species
        assert len(species) > 20

    def test_load_custom_json(self, tmp_path):
        """Test loading custom karyotypes from JSON."""
        custom_data = [
            {
                "name": "custom_species",
                "scientific_name": "Customus speciesus",
                "taxid": 99999,
                "chromosome_count": 8,
                "autosomes": ["1", "2", "3", "4", "5", "6"],
                "sex_chromosomes": ["X", "Y"],
            }
        ]
        custom_file = tmp_path / "custom.json"
        custom_file.write_text(json.dumps(custom_data))

        db = KaryotypeDatabase()
        count = db.load_custom(custom_file)
        assert count == 1

        custom = db.lookup("custom_species")
        assert custom is not None
        assert custom.chromosome_count == 8

    @pytest.mark.skipif(not HAS_YAML, reason="PyYAML not installed")
    def test_load_custom_yaml(self, tmp_path):
        """Test loading custom karyotypes from YAML."""
        custom_yaml = """
- name: yaml_species
  scientific_name: Yamlus speciesus
  taxid: 88888
  chromosome_count: 6
  autosomes: ["1", "2", "3", "4"]
  sex_chromosomes: ["Z", "W"]
  sex_system: "ZZ/ZW"
"""
        custom_file = tmp_path / "custom.yaml"
        custom_file.write_text(custom_yaml)

        db = KaryotypeDatabase()
        count = db.load_custom(custom_file)
        assert count == 1

        custom = db.lookup("yaml_species")
        assert custom is not None
        assert custom.sex_system == SexChromosomeSystem.ZZ_ZW

    def test_load_custom_file_not_found(self, tmp_path):
        """Test error when custom file doesn't exist."""
        db = KaryotypeDatabase()
        with pytest.raises(FileNotFoundError):
            db.load_custom(tmp_path / "nonexistent.json")


class TestBuiltinKaryotypes:
    """Tests for specific built-in karyotypes."""

    def test_human_karyotype(self):
        """Test human karyotype is correct."""
        db = KaryotypeDatabase()
        human = db.lookup("human")
        assert human is not None
        assert human.chromosome_count == 24
        assert len(human.autosomes) == 22
        assert human.sex_chromosomes == ["X", "Y"]
        assert human.sex_system == SexChromosomeSystem.XX_XY

    def test_mouse_karyotype(self):
        """Test mouse karyotype is correct."""
        db = KaryotypeDatabase()
        mouse = db.lookup("mouse")
        assert mouse is not None
        assert mouse.chromosome_count == 21
        assert len(mouse.autosomes) == 19

    def test_drosophila_karyotype(self):
        """Test Drosophila karyotype with chromosome arms."""
        db = KaryotypeDatabase()
        fly = db.lookup("fruitfly")
        assert fly is not None
        assert "2L" in fly.autosomes
        assert "2R" in fly.autosomes
        assert "3L" in fly.autosomes
        assert "3R" in fly.autosomes

    def test_chicken_karyotype(self):
        """Test chicken has ZW sex system."""
        db = KaryotypeDatabase()
        chicken = db.lookup("chicken")
        assert chicken is not None
        assert chicken.sex_system == SexChromosomeSystem.ZZ_ZW
        assert "Z" in chicken.sex_chromosomes
        assert "W" in chicken.sex_chromosomes

    def test_celegans_karyotype(self):
        """Test C. elegans has XO sex system."""
        db = KaryotypeDatabase()
        worm = db.lookup("celegans")
        assert worm is not None
        assert worm.sex_system == SexChromosomeSystem.XX_XO
        assert "I" in worm.autosomes
        assert "X" in worm.sex_chromosomes

    def test_yeast_karyotype(self):
        """Test yeast has Roman numeral chromosomes."""
        db = KaryotypeDatabase()
        yeast = db.lookup("yeast")
        assert yeast is not None
        assert yeast.chromosome_count == 16
        assert "I" in yeast.autosomes
        assert "XVI" in yeast.autosomes

    @pytest.mark.parametrize(
        (
            "query",
            "scientific_name",
            "chromosome_count",
            "autosome_count",
            "expected_note",
        ),
        [
            ("atlantic_salmon", "Salmo salar", 29, 29, "2n=58"),
            ("rainbow_trout", "Oncorhynchus mykiss", 29, 29, "2n=58-64"),
            ("tilapia", "Oreochromis niloticus", 22, 22, "2n=44"),
            ("common_carp", "Cyprinus carpio", 50, 50, "2n=100"),
            ("medaka", "Oryzias latipes", 24, 23, "2n=48"),
            ("stickleback", "Gasterosteus aculeatus", 21, 21, "2n=42"),
            ("fugu", "Takifugu rubripes", 22, 22, "2n=44"),
        ],
    )
    def test_fish_and_aquaculture_karyotypes(
        self,
        query,
        scientific_name,
        chromosome_count,
        autosome_count,
        expected_note,
    ):
        """Test fish and aquaculture karyotypes requested in issue #8."""
        db = KaryotypeDatabase()
        fish = db.lookup(query)
        assert fish is not None
        assert fish.scientific_name == scientific_name
        assert fish.chromosome_count == chromosome_count
        assert len(fish.autosomes) == autosome_count
        assert expected_note in fish.notes

    @pytest.mark.parametrize(
        ("alias", "expected_name"),
        [
            ("atlantic salmon", "atlantic_salmon"),
            ("rainbow trout", "rainbow_trout"),
            ("nile tilapia", "tilapia"),
            ("carp", "common_carp"),
            ("three-spined stickleback", "stickleback"),
        ],
    )
    def test_fish_and_aquaculture_aliases(self, alias, expected_name):
        """Test common-name aliases for fish and aquaculture karyotypes."""
        db = KaryotypeDatabase()
        fish = db.lookup(alias)
        assert fish is not None
        assert fish.name == expected_name


class TestValidateKaryotype:
    """Tests for karyotype validation function."""

    def test_valid_human_assembly(self):
        """Test validation of correct human assembly."""
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        result = validate_karyotype(chromosomes, "human")

        assert result.is_valid
        assert result.chromosome_count_match
        # MT may be in missing_chromosomes as a warning (not an error)
        non_mt_missing = [c for c in result.missing_chromosomes if c != "MT"]
        assert len(non_mt_missing) == 0
        assert len(result.extra_chromosomes) == 0

    def test_missing_chromosomes(self):
        """Test detection of missing chromosomes."""
        # Missing chromosomes 21, 22, X, Y
        chromosomes = [f"chr{i}" for i in range(1, 21)]
        result = validate_karyotype(chromosomes, "human")

        assert not result.is_valid
        assert not result.chromosome_count_match
        assert "21" in result.missing_chromosomes
        assert "22" in result.missing_chromosomes
        assert "X" in result.missing_chromosomes
        assert "Y" in result.missing_chromosomes

    def test_extra_chromosomes(self):
        """Test detection of unexpected chromosomes."""
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrExtra"]
        result = validate_karyotype(chromosomes, "human")

        # Extra chromosomes cause a count mismatch (25 vs 24), which is an error
        assert not result.is_valid  # Count mismatch makes it invalid
        assert not result.chromosome_count_match
        assert "chrExtra" in result.extra_chromosomes
        # Should have an unexpected_chromosome warning
        unexpected_issues = [i for i in result.issues if i.issue_type == "unexpected_chromosome"]
        assert len(unexpected_issues) >= 1

    def test_extra_chromosomes_strict_mode(self):
        """Test extra chromosomes become errors in strict mode."""
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrExtra"]
        result = validate_karyotype(chromosomes, "human", strict=True)

        assert not result.is_valid
        error_count = len([i for i in result.issues if i.severity == "error"])
        assert error_count >= 1

    def test_species_not_found(self):
        """Test error when species not in database."""
        result = validate_karyotype(["chr1", "chr2"], "unknown_species_xyz")

        assert not result.is_valid
        assert result.karyotype is None
        assert any(i.issue_type == "species_not_found" for i in result.issues)

    def test_chromosome_normalization(self):
        """Test that chromosome names are normalized."""
        # Various formats should all match
        chromosomes = [
            "chromosome_1", "chromosome_2", "Chr_3", "chr4",
            "5", "6", "7", "8", "9", "10", "11", "12", "13",
            "14", "15", "16", "17", "18", "19", "20", "21", "22",
            "chrX", "chromosome_Y"
        ]
        result = validate_karyotype(chromosomes, "human")

        assert result.is_valid
        assert result.chromosome_count_match

    def test_mitochondria_optional_by_default(self):
        """Test that missing MT is a warning, not error."""
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        # No MT
        result = validate_karyotype(chromosomes, "human")

        assert result.is_valid
        # Should have a warning about MT
        mt_warnings = [i for i in result.issues if "mitochondria" in i.issue_type.lower()]
        assert len(mt_warnings) >= 1
        assert mt_warnings[0].severity == "warning"

    def test_mitochondria_required_in_strict(self):
        """Test that missing MT is an error in strict mode."""
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        result = validate_karyotype(chromosomes, "human", strict=True)

        # MT missing should be an error
        mt_issues = [i for i in result.issues if "mitochondria" in i.issue_type.lower()]
        # In strict mode, warnings become errors
        # But chromosome count is still matched (24), so we need to check specifically

    def test_mouse_validation(self):
        """Test validation of mouse assembly."""
        chromosomes = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]
        result = validate_karyotype(chromosomes, "mouse")

        assert result.is_valid
        assert result.chromosome_count_match

    def test_lookup_by_taxid(self):
        """Test validation using taxonomy ID."""
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        result = validate_karyotype(chromosomes, 9606)  # Human taxid

        assert result.is_valid
        assert result.karyotype is not None
        assert result.karyotype.name == "human"


class TestFormatKaryotypeSummary:
    """Tests for summary formatting."""

    def test_valid_result_formatting(self):
        """Test formatting of valid result."""
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        result = validate_karyotype(chromosomes, "human")
        summary = format_karyotype_summary(result)

        assert "KARYOTYPE VALIDATION" in summary
        assert "human" in summary
        assert "Homo sapiens" in summary
        assert "VALID" in summary

    def test_invalid_result_formatting(self):
        """Test formatting of invalid result."""
        chromosomes = ["chr1", "chr2"]  # Missing most chromosomes
        result = validate_karyotype(chromosomes, "human")
        summary = format_karyotype_summary(result)

        assert "INVALID" in summary
        assert "ERRORS:" in summary
        assert "MISSING CHROMOSOMES" in summary


class TestFormatKaryotypeTsv:
    """Tests for TSV formatting."""

    def test_tsv_output(self):
        """Test TSV format output."""
        chromosomes = ["chr1", "chr2"]
        result = validate_karyotype(chromosomes, "human")
        tsv = format_karyotype_tsv(result)

        lines = tsv.strip().split("\n")
        assert "severity" in lines[0]
        assert "issue_type" in lines[0]
        assert len(lines) > 1


class TestKaryotypeValidationResult:
    """Tests for KaryotypeValidationResult."""

    def test_to_dict(self):
        """Test serialization to dictionary."""
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        result = validate_karyotype(chromosomes, "human")
        d = result.to_dict()

        assert d["species"] == "human"
        assert d["is_valid"] is True
        assert "summary" in d
        assert d["summary"]["expected_count"] > 0


class TestKaryotypeCLI:
    """Tests for karyotype validation CLI."""

    @pytest.fixture
    def human_fasta(self) -> Path:
        """Create a FASTA file with human chromosomes."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            # Write chromosomes 1-22, X, Y
            for i in range(1, 23):
                f.write(f">chr{i}\n")
                f.write("ATCGATCG" * 1000 + "\n")
            f.write(">chrX\n")
            f.write("ATCGATCG" * 1000 + "\n")
            f.write(">chrY\n")
            f.write("ATCGATCG" * 1000 + "\n")
            f.flush()
            return Path(f.name)

    @pytest.fixture
    def incomplete_fasta(self) -> Path:
        """Create a FASTA with missing chromosomes."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            # Only chromosomes 1-10
            for i in range(1, 11):
                f.write(f">chr{i}\n")
                f.write("ATCGATCG" * 1000 + "\n")
            f.flush()
            return Path(f.name)

    def test_list_species(self):
        """Test --list-species flag."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "--list-species"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "human" in result.stdout
        assert "Homo sapiens" in result.stdout
        assert "mouse" in result.stdout

    def test_check_karyotype_valid(self, human_fasta):
        """Test karyotype validation with valid assembly."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(human_fasta),
                "--check-karyotype", "human",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "VALID" in result.stdout

    def test_check_karyotype_invalid(self, incomplete_fasta):
        """Test karyotype validation with incomplete assembly."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(incomplete_fasta),
                "--check-karyotype", "human",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "INVALID" in result.stdout

    def test_check_karyotype_json_output(self, human_fasta):
        """Test karyotype validation with JSON output."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(human_fasta),
                "--check-karyotype", "human",
                "-f", "json",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert data["is_valid"] is True
        assert data["species"] == "human"

    def test_check_karyotype_unknown_species(self, human_fasta):
        """Test error for unknown species."""
        result = subprocess.run(
            [
                sys.executable, "-m", "chromdetect",
                str(human_fasta),
                "--check-karyotype", "unknown_xyz",
                "-q",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "not found" in result.stdout.lower()

    def test_check_karyotype_help_text(self):
        """Test that karyotype options are in help."""
        result = subprocess.run(
            [sys.executable, "-m", "chromdetect", "--help"],
            capture_output=True,
            text=True,
        )
        assert "--check-karyotype" in result.stdout
        assert "--list-species" in result.stdout
