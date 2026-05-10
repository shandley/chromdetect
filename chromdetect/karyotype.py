"""
Karyotype Database and Validation for ChromDetect.

This module provides a database of expected karyotypes for common species
and validates genome assemblies against these expectations.

The karyotype database includes:
- Expected chromosome count
- Expected chromosome names/IDs
- Sex chromosome system information
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

try:
    import yaml  # type: ignore[import-untyped]
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


class SexChromosomeSystem(Enum):
    """Sex chromosome determination systems."""

    XX_XY = "XX/XY"  # Female XX, Male XY (mammals, some insects)
    ZZ_ZW = "ZZ/ZW"  # Male ZZ, Female ZW (birds, some reptiles, fish)
    XX_XO = "XX/XO"  # Female XX, Male XO (some insects)
    HAPLODIPLOID = "haplodiploid"  # Females diploid, males haploid (bees, ants)
    NONE = "none"  # No sex chromosomes (some plants, fungi)
    UNKNOWN = "unknown"


@dataclass
class Karyotype:
    """Expected karyotype information for a species.

    Attributes:
        name: Common name of the species
        scientific_name: Scientific (Latin) name
        taxid: NCBI Taxonomy ID
        chromosome_count: Total expected chromosome count (haploid n)
        autosomes: List of autosome IDs (e.g., ["1", "2", ..., "22"])
        sex_chromosomes: List of sex chromosome IDs (e.g., ["X", "Y"])
        sex_system: Sex chromosome determination system
        mitochondria: Whether mitochondrial genome is expected
        aliases: Alternative names for lookup
        notes: Additional notes about the karyotype
    """

    name: str
    scientific_name: str
    taxid: int | None
    chromosome_count: int
    autosomes: list[str]
    sex_chromosomes: list[str] = field(default_factory=list)
    sex_system: SexChromosomeSystem = SexChromosomeSystem.UNKNOWN
    mitochondria: bool = True
    aliases: list[str] = field(default_factory=list)
    notes: str = ""

    @property
    def all_chromosomes(self) -> list[str]:
        """Get all expected chromosome IDs."""
        return self.autosomes + self.sex_chromosomes

    @property
    def total_with_mt(self) -> int:
        """Get total count including mitochondria."""
        return self.chromosome_count + (1 if self.mitochondria else 0)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "scientific_name": self.scientific_name,
            "taxid": self.taxid,
            "chromosome_count": self.chromosome_count,
            "autosomes": self.autosomes,
            "sex_chromosomes": self.sex_chromosomes,
            "sex_system": self.sex_system.value,
            "mitochondria": self.mitochondria,
            "aliases": self.aliases,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict) -> Karyotype:
        """Create Karyotype from dictionary."""
        sex_system_str = data.get("sex_system", "unknown")
        try:
            sex_system = SexChromosomeSystem(sex_system_str)
        except ValueError:
            sex_system = SexChromosomeSystem.UNKNOWN

        return cls(
            name=data["name"],
            scientific_name=data.get("scientific_name", ""),
            taxid=data.get("taxid"),
            chromosome_count=data["chromosome_count"],
            autosomes=data.get("autosomes", []),
            sex_chromosomes=data.get("sex_chromosomes", []),
            sex_system=sex_system,
            mitochondria=data.get("mitochondria", True),
            aliases=data.get("aliases", []),
            notes=data.get("notes", ""),
        )


# Built-in karyotype database for common model organisms and species
KARYOTYPE_DATABASE: list[dict[str, Any]] = [
    # Mammals
    {
        "name": "human",
        "scientific_name": "Homo sapiens",
        "taxid": 9606,
        "chromosome_count": 24,
        "autosomes": [str(i) for i in range(1, 23)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["homo sapiens", "h. sapiens", "hsapiens"],
    },
    {
        "name": "mouse",
        "scientific_name": "Mus musculus",
        "taxid": 10090,
        "chromosome_count": 21,
        "autosomes": [str(i) for i in range(1, 20)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["house mouse", "m. musculus", "mmusculus"],
    },
    {
        "name": "rat",
        "scientific_name": "Rattus norvegicus",
        "taxid": 10116,
        "chromosome_count": 22,
        "autosomes": [str(i) for i in range(1, 21)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["norway rat", "brown rat", "r. norvegicus", "rnorvegicus"],
    },
    {
        "name": "dog",
        "scientific_name": "Canis lupus familiaris",
        "taxid": 9615,
        "chromosome_count": 40,
        "autosomes": [str(i) for i in range(1, 39)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["domestic dog", "canis familiaris", "cfamiliaris"],
    },
    {
        "name": "cat",
        "scientific_name": "Felis catus",
        "taxid": 9685,
        "chromosome_count": 20,
        "autosomes": ["A1", "A2", "A3", "B1", "B2", "B3", "B4", "C1", "C2", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "F1", "F2"],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["domestic cat", "f. catus", "fcatus"],
        "notes": "Cat chromosomes use letter-number nomenclature (A1-F2)",
    },
    {
        "name": "pig",
        "scientific_name": "Sus scrofa",
        "taxid": 9823,
        "chromosome_count": 20,
        "autosomes": [str(i) for i in range(1, 19)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["domestic pig", "wild boar", "s. scrofa", "sscrofa"],
    },
    {
        "name": "cow",
        "scientific_name": "Bos taurus",
        "taxid": 9913,
        "chromosome_count": 31,
        "autosomes": [str(i) for i in range(1, 30)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["cattle", "domestic cattle", "b. taurus", "btaurus"],
    },
    {
        "name": "horse",
        "scientific_name": "Equus caballus",
        "taxid": 9796,
        "chromosome_count": 33,
        "autosomes": [str(i) for i in range(1, 32)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["domestic horse", "e. caballus", "ecaballus"],
    },
    {
        "name": "sheep",
        "scientific_name": "Ovis aries",
        "taxid": 9940,
        "chromosome_count": 28,
        "autosomes": [str(i) for i in range(1, 27)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["domestic sheep", "o. aries", "oaries"],
    },
    {
        "name": "goat",
        "scientific_name": "Capra hircus",
        "taxid": 9925,
        "chromosome_count": 31,
        "autosomes": [str(i) for i in range(1, 30)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["domestic goat", "c. hircus", "chircus"],
    },
    # Birds
    {
        "name": "chicken",
        "scientific_name": "Gallus gallus",
        "taxid": 9031,
        "chromosome_count": 40,
        "autosomes": [str(i) for i in range(1, 39)],
        "sex_chromosomes": ["Z", "W"],
        "sex_system": "ZZ/ZW",
        "aliases": ["red junglefowl", "g. gallus", "ggallus"],
        "notes": "Many microchromosomes; count may vary by assembly",
    },
    {
        "name": "zebrafinch",
        "scientific_name": "Taeniopygia guttata",
        "taxid": 59729,
        "chromosome_count": 40,
        "autosomes": [str(i) for i in range(1, 39)],
        "sex_chromosomes": ["Z", "W"],
        "sex_system": "ZZ/ZW",
        "aliases": ["zebra finch", "t. guttata", "tguttata"],
    },
    # Fish
    {
        "name": "zebrafish",
        "scientific_name": "Danio rerio",
        "taxid": 7955,
        "chromosome_count": 25,
        "autosomes": [str(i) for i in range(1, 26)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["d. rerio", "drerio"],
        "notes": "No differentiated sex chromosomes",
    },
    {
        "name": "medaka",
        "scientific_name": "Oryzias latipes",
        "taxid": 8090,
        "chromosome_count": 24,
        "autosomes": [str(i) for i in range(1, 24)],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["japanese medaka", "japanese rice fish", "o. latipes", "olatipes"],
    },
    {
        "name": "fugu",
        "scientific_name": "Takifugu rubripes",
        "taxid": 31033,
        "chromosome_count": 22,
        "autosomes": [str(i) for i in range(1, 23)],
        "sex_chromosomes": [],
        "sex_system": "unknown",
        "aliases": ["japanese pufferfish", "t. rubripes", "trubripes"],
    },
    # Amphibians
    {
        "name": "xenopus",
        "scientific_name": "Xenopus tropicalis",
        "taxid": 8364,
        "chromosome_count": 10,
        "autosomes": [str(i) for i in range(1, 10)],
        "sex_chromosomes": ["W", "Z"],
        "sex_system": "ZZ/ZW",
        "aliases": ["western clawed frog", "x. tropicalis", "xtropicalis"],
    },
    # Insects
    {
        "name": "fruitfly",
        "scientific_name": "Drosophila melanogaster",
        "taxid": 7227,
        "chromosome_count": 5,
        "autosomes": ["2L", "2R", "3L", "3R", "4"],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["fruit fly", "drosophila", "d. melanogaster", "dmelanogaster"],
        "notes": "Chromosome arms often listed separately (2L, 2R, 3L, 3R)",
    },
    {
        "name": "mosquito",
        "scientific_name": "Anopheles gambiae",
        "taxid": 7165,
        "chromosome_count": 5,
        "autosomes": ["2L", "2R", "3L", "3R"],
        "sex_chromosomes": ["X", "Y"],
        "sex_system": "XX/XY",
        "aliases": ["african malaria mosquito", "a. gambiae", "agambiae"],
    },
    {
        "name": "honeybee",
        "scientific_name": "Apis mellifera",
        "taxid": 7460,
        "chromosome_count": 16,
        "autosomes": [str(i) for i in range(1, 17)],
        "sex_chromosomes": [],
        "sex_system": "haplodiploid",
        "aliases": ["western honeybee", "a. mellifera", "amellifera"],
        "notes": "Haplodiploid: males haploid (16), females diploid (32)",
    },
    # Nematodes
    {
        "name": "celegans",
        "scientific_name": "Caenorhabditis elegans",
        "taxid": 6239,
        "chromosome_count": 6,
        "autosomes": ["I", "II", "III", "IV", "V"],
        "sex_chromosomes": ["X"],
        "sex_system": "XX/XO",
        "aliases": ["c. elegans", "caenorhabditis elegans", "roundworm"],
        "notes": "Hermaphrodites XX, males XO",
    },
    # Plants
    {
        "name": "arabidopsis",
        "scientific_name": "Arabidopsis thaliana",
        "taxid": 3702,
        "chromosome_count": 5,
        "autosomes": ["1", "2", "3", "4", "5"],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["thale cress", "a. thaliana", "athaliana"],
    },
    {
        "name": "rice",
        "scientific_name": "Oryza sativa",
        "taxid": 4530,
        "chromosome_count": 12,
        "autosomes": [str(i) for i in range(1, 13)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["asian rice", "o. sativa", "osativa"],
    },
    {
        "name": "maize",
        "scientific_name": "Zea mays",
        "taxid": 4577,
        "chromosome_count": 10,
        "autosomes": [str(i) for i in range(1, 11)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["corn", "z. mays", "zmays"],
    },
    {
        "name": "wheat",
        "scientific_name": "Triticum aestivum",
        "taxid": 4565,
        "chromosome_count": 21,
        "autosomes": [f"{i}{g}" for g in ["A", "B", "D"] for i in range(1, 8)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["bread wheat", "common wheat", "t. aestivum", "taestivum"],
        "notes": "Hexaploid with A, B, D subgenomes",
    },
    {
        "name": "tomato",
        "scientific_name": "Solanum lycopersicum",
        "taxid": 4081,
        "chromosome_count": 12,
        "autosomes": [str(i) for i in range(1, 13)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["s. lycopersicum", "slycopersicum"],
    },
    {
        "name": "soybean",
        "scientific_name": "Glycine max",
        "taxid": 3847,
        "chromosome_count": 20,
        "autosomes": [str(i) for i in range(1, 21)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["soya bean", "g. max", "gmax"],
    },
    {
        "name": "potato",
        "scientific_name": "Solanum tuberosum",
        "taxid": 4113,
        "chromosome_count": 24,
        "autosomes": [str(i) for i in range(1, 25)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["s. tuberosum", "stuberosum"],
        "notes": "Tetraploid crop species, commonly reported as 2n=48",
    },
    {
        "name": "cotton",
        "scientific_name": "Gossypium hirsutum",
        "taxid": 3635,
        "chromosome_count": 26,
        "autosomes": [str(i) for i in range(1, 27)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["upland cotton", "g. hirsutum", "ghirsutum"],
        "notes": "Allotetraploid crop species, commonly reported as 2n=52",
    },
    {
        "name": "sunflower",
        "scientific_name": "Helianthus annuus",
        "taxid": 4232,
        "chromosome_count": 17,
        "autosomes": [str(i) for i in range(1, 18)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["common sunflower", "h. annuus", "hannuus"],
        "notes": "Commonly reported as 2n=34",
    },
    {
        "name": "barley",
        "scientific_name": "Hordeum vulgare",
        "taxid": 4513,
        "chromosome_count": 7,
        "autosomes": [str(i) for i in range(1, 8)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["h. vulgare", "hvulgare"],
        "notes": "Commonly reported as 2n=14",
    },
    {
        "name": "grape",
        "scientific_name": "Vitis vinifera",
        "taxid": 29760,
        "chromosome_count": 19,
        "autosomes": [str(i) for i in range(1, 20)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["grapevine", "v. vinifera", "vvinifera"],
        "notes": "Commonly reported as 2n=38",
    },
    {
        "name": "cannabis",
        "scientific_name": "Cannabis sativa",
        "taxid": 3483,
        "chromosome_count": 10,
        "autosomes": [str(i) for i in range(1, 11)],
        "sex_chromosomes": [],
        "sex_system": "unknown",
        "aliases": ["hemp", "c. sativa", "csativa"],
        "notes": "Dioecious species, commonly reported as 2n=20",
    },
    {
        "name": "strawberry",
        "scientific_name": "Fragaria × ananassa",
        "taxid": 3747,
        "chromosome_count": 28,
        "autosomes": [str(i) for i in range(1, 29)],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["garden strawberry", "fragaria ananassa", "f. ananassa", "fananassa"],
        "notes": "Octoploid crop species, commonly reported as 2n=56",
    },
    # Fungi
    {
        "name": "yeast",
        "scientific_name": "Saccharomyces cerevisiae",
        "taxid": 4932,
        "chromosome_count": 16,
        "autosomes": ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"],
        "sex_chromosomes": [],
        "sex_system": "none",
        "mitochondria": True,
        "aliases": ["budding yeast", "baker's yeast", "s. cerevisiae", "scerevisiae"],
    },
    {
        "name": "fission_yeast",
        "scientific_name": "Schizosaccharomyces pombe",
        "taxid": 4896,
        "chromosome_count": 3,
        "autosomes": ["I", "II", "III"],
        "sex_chromosomes": [],
        "sex_system": "none",
        "aliases": ["s. pombe", "spombe"],
    },
]


class KaryotypeDatabase:
    """Database of karyotype information for species lookup."""

    def __init__(self) -> None:
        """Initialize with built-in database."""
        self._karyotypes: dict[str, Karyotype] = {}
        self._by_taxid: dict[int, Karyotype] = {}
        self._load_builtin()

    def _load_builtin(self) -> None:
        """Load the built-in karyotype database."""
        for entry in KARYOTYPE_DATABASE:
            karyotype = Karyotype.from_dict(entry)
            self._add_karyotype(karyotype)

    def _add_karyotype(self, karyotype: Karyotype) -> None:
        """Add a karyotype to the database with all lookup keys."""
        # Primary name (lowercase)
        self._karyotypes[karyotype.name.lower()] = karyotype

        # Scientific name (lowercase)
        if karyotype.scientific_name:
            self._karyotypes[karyotype.scientific_name.lower()] = karyotype

        # Aliases
        for alias in karyotype.aliases:
            self._karyotypes[alias.lower()] = karyotype

        # Taxonomy ID
        if karyotype.taxid:
            self._by_taxid[karyotype.taxid] = karyotype

    def load_custom(self, file_path: Path | str) -> int:
        """
        Load additional karyotypes from a JSON or YAML file.

        Args:
            file_path: Path to custom karyotype file

        Returns:
            Number of karyotypes loaded

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"Karyotype file not found: {file_path}")

        content = file_path.read_text()
        ext = file_path.suffix.lower()

        if ext in (".yaml", ".yml"):
            if not HAS_YAML:
                raise ValueError(
                    "YAML support requires PyYAML. Install with: pip install pyyaml"
                )
            try:
                data = yaml.safe_load(content)
            except yaml.YAMLError as e:
                raise ValueError(f"Invalid YAML: {e}") from e
        else:
            try:
                data = json.loads(content)
            except json.JSONDecodeError as e:
                raise ValueError(f"Invalid JSON: {e}") from e

        if not isinstance(data, list):
            if "karyotypes" in data:
                data = data["karyotypes"]
            else:
                raise ValueError("Expected list of karyotypes or object with 'karyotypes' key")

        count = 0
        for entry in data:
            karyotype = Karyotype.from_dict(entry)
            self._add_karyotype(karyotype)
            count += 1

        return count

    def lookup(self, query: str | int) -> Karyotype | None:
        """
        Look up a karyotype by name, scientific name, alias, or taxid.

        Args:
            query: Species name (string) or NCBI taxonomy ID (int)

        Returns:
            Karyotype if found, None otherwise
        """
        if isinstance(query, int):
            return self._by_taxid.get(query)
        return self._karyotypes.get(query.lower())

    def list_species(self) -> list[str]:
        """Get list of all available species names."""
        seen = set()
        names = []
        for karyotype in self._karyotypes.values():
            if karyotype.name not in seen:
                seen.add(karyotype.name)
                names.append(karyotype.name)
        return sorted(names)

    def __len__(self) -> int:
        """Return number of unique species in database."""
        return len({id(k) for k in self._karyotypes.values()})


@dataclass
class KaryotypeValidationIssue:
    """An issue found during karyotype validation."""

    severity: str  # "error", "warning", "info"
    issue_type: str
    message: str
    expected: str | int | None = None
    found: str | int | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "severity": self.severity,
            "issue_type": self.issue_type,
            "message": self.message,
            "expected": self.expected,
            "found": self.found,
        }


@dataclass
class KaryotypeValidationResult:
    """Result of validating an assembly against an expected karyotype.

    Attributes:
        species: Species name used for lookup
        karyotype: The expected karyotype (None if species not found)
        is_valid: True if assembly matches expected karyotype
        issues: List of validation issues
        found_chromosomes: Chromosomes found in assembly
        expected_chromosomes: Expected chromosomes from karyotype
        missing_chromosomes: Expected but not found
        extra_chromosomes: Found but not expected
        chromosome_count_match: Whether total counts match
    """

    species: str
    karyotype: Karyotype | None
    is_valid: bool
    issues: list[KaryotypeValidationIssue] = field(default_factory=list)
    found_chromosomes: list[str] = field(default_factory=list)
    expected_chromosomes: list[str] = field(default_factory=list)
    missing_chromosomes: list[str] = field(default_factory=list)
    extra_chromosomes: list[str] = field(default_factory=list)
    chromosome_count_match: bool = False

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "species": self.species,
            "karyotype": self.karyotype.to_dict() if self.karyotype else None,
            "is_valid": self.is_valid,
            "summary": {
                "expected_count": len(self.expected_chromosomes),
                "found_count": len(self.found_chromosomes),
                "missing_count": len(self.missing_chromosomes),
                "extra_count": len(self.extra_chromosomes),
                "chromosome_count_match": self.chromosome_count_match,
            },
            "issues": [i.to_dict() for i in self.issues],
            "found_chromosomes": self.found_chromosomes,
            "expected_chromosomes": self.expected_chromosomes,
            "missing_chromosomes": self.missing_chromosomes,
            "extra_chromosomes": self.extra_chromosomes,
        }


def _normalize_chromosome_id(chrom_id: str) -> str:
    """Normalize chromosome ID for comparison.

    Handles variations like:
    - chr1 -> 1
    - chromosome_1 -> 1
    - Chr_X -> X
    - chrMT -> MT
    """
    normalized = chrom_id.lower()

    # Remove common prefixes
    for prefix in ["chromosome_", "chromosome", "chr_", "chr"]:
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix):]
            break

    return normalized.upper()


def validate_karyotype(
    assembly_chromosomes: list[str],
    species: str | int,
    database: KaryotypeDatabase | None = None,
    strict: bool = False,
) -> KaryotypeValidationResult:
    """
    Validate assembly chromosomes against expected karyotype.

    Args:
        assembly_chromosomes: List of chromosome names/IDs from assembly
        species: Species name or NCBI taxonomy ID
        database: KaryotypeDatabase instance (uses default if None)
        strict: If True, treat warnings as errors

    Returns:
        KaryotypeValidationResult with findings
    """
    if database is None:
        database = KaryotypeDatabase()

    karyotype = database.lookup(species)
    issues: list[KaryotypeValidationIssue] = []

    if karyotype is None:
        return KaryotypeValidationResult(
            species=str(species),
            karyotype=None,
            is_valid=False,
            issues=[KaryotypeValidationIssue(
                severity="error",
                issue_type="species_not_found",
                message=f"Species '{species}' not found in karyotype database",
            )],
        )

    # Normalize chromosome IDs for comparison
    found_normalized = {_normalize_chromosome_id(c): c for c in assembly_chromosomes}
    expected_normalized = {_normalize_chromosome_id(c): c for c in karyotype.all_chromosomes}

    # Also add MT/M normalization
    if karyotype.mitochondria:
        expected_normalized["MT"] = "MT"
        expected_normalized["M"] = "MT"

    # Find matches, missing, and extra
    found_set = set(found_normalized.keys())
    expected_set = set(expected_normalized.keys())

    # Handle MT/M as equivalent
    if "MT" in found_set or "M" in found_set:
        found_set.discard("M")
        if "M" in found_normalized and "MT" not in found_normalized:
            found_normalized["MT"] = found_normalized["M"]
        found_set.add("MT")

    missing = expected_set - found_set
    extra = found_set - expected_set

    # Remove MT from missing if not strict (it's optional)
    if not strict and "MT" in missing:
        missing.discard("MT")

    # Check chromosome count
    expected_count = karyotype.chromosome_count
    found_count = len([c for c in found_set if c not in ("MT", "M")])
    count_match = found_count == expected_count

    # Build issues list
    if not count_match:
        issues.append(KaryotypeValidationIssue(
            severity="error",
            issue_type="chromosome_count_mismatch",
            message=f"Expected {expected_count} chromosomes, found {found_count}",
            expected=expected_count,
            found=found_count,
        ))

    for chrom in sorted(missing):
        if chrom in ("MT", "M"):
            severity = "warning" if not strict else "error"
            issues.append(KaryotypeValidationIssue(
                severity=severity,
                issue_type="missing_mitochondria",
                message="Mitochondrial genome (MT) not found in assembly",
                expected="MT",
            ))
        else:
            issues.append(KaryotypeValidationIssue(
                severity="error",
                issue_type="missing_chromosome",
                message=f"Expected chromosome '{chrom}' not found in assembly",
                expected=chrom,
            ))

    for chrom in sorted(extra):
        original_name = found_normalized.get(chrom, chrom)
        issues.append(KaryotypeValidationIssue(
            severity="warning" if not strict else "error",
            issue_type="unexpected_chromosome",
            message=f"Unexpected chromosome '{original_name}' (normalized: {chrom}) not in expected karyotype",
            found=original_name,
        ))

    # Determine validity
    error_count = len([i for i in issues if i.severity == "error"])
    is_valid = error_count == 0

    return KaryotypeValidationResult(
        species=str(species),
        karyotype=karyotype,
        is_valid=is_valid,
        issues=issues,
        found_chromosomes=assembly_chromosomes,
        expected_chromosomes=karyotype.all_chromosomes + (["MT"] if karyotype.mitochondria else []),
        missing_chromosomes=[expected_normalized.get(c, c) for c in sorted(missing)],
        extra_chromosomes=[found_normalized.get(c, c) for c in sorted(extra)],
        chromosome_count_match=count_match,
    )


def format_karyotype_summary(result: KaryotypeValidationResult) -> str:
    """
    Format karyotype validation result as human-readable summary.

    Args:
        result: KaryotypeValidationResult

    Returns:
        Formatted string for display
    """
    lines = [
        "=" * 60,
        "KARYOTYPE VALIDATION",
        "=" * 60,
        "",
        f"Species query:   {result.species}",
    ]

    if result.karyotype:
        lines.extend([
            f"Matched species: {result.karyotype.name}",
            f"Scientific name: {result.karyotype.scientific_name}",
            f"Expected chromosomes: {result.karyotype.chromosome_count}",
            f"Sex system:      {result.karyotype.sex_system.value}",
        ])
    else:
        lines.append("Matched species: NOT FOUND")

    lines.extend([
        "",
        "-" * 60,
        "SUMMARY",
        "-" * 60,
        f"Expected chromosomes: {len(result.expected_chromosomes)}",
        f"Found chromosomes:    {len(result.found_chromosomes)}",
        f"Missing:              {len(result.missing_chromosomes)}",
        f"Unexpected:           {len(result.extra_chromosomes)}",
        "",
    ])

    if result.is_valid:
        lines.append("STATUS: VALID")
    else:
        error_count = len([i for i in result.issues if i.severity == "error"])
        lines.append(f"STATUS: INVALID - {error_count} error(s) found")

    warning_count = len([i for i in result.issues if i.severity == "warning"])
    if warning_count > 0:
        lines.append(f"         {warning_count} warning(s)")

    if result.issues:
        lines.extend([
            "",
            "-" * 60,
            "ISSUES",
            "-" * 60,
        ])

        errors = [i for i in result.issues if i.severity == "error"]
        warnings = [i for i in result.issues if i.severity == "warning"]

        if errors:
            lines.append("")
            lines.append("ERRORS:")
            for issue in errors:
                lines.append(f"  [{issue.issue_type}] {issue.message}")

        if warnings:
            lines.append("")
            lines.append("WARNINGS:")
            for issue in warnings:
                lines.append(f"  [{issue.issue_type}] {issue.message}")

    if result.missing_chromosomes:
        lines.extend([
            "",
            "-" * 60,
            f"MISSING CHROMOSOMES ({len(result.missing_chromosomes)}):",
            "-" * 60,
            "  " + ", ".join(result.missing_chromosomes),
        ])

    if result.extra_chromosomes:
        lines.extend([
            "",
            "-" * 60,
            f"UNEXPECTED CHROMOSOMES ({len(result.extra_chromosomes)}):",
            "-" * 60,
            "  " + ", ".join(result.extra_chromosomes),
        ])

    return "\n".join(lines)


def format_karyotype_tsv(result: KaryotypeValidationResult) -> str:
    """
    Format karyotype validation issues as TSV.

    Args:
        result: KaryotypeValidationResult

    Returns:
        TSV-formatted string
    """
    lines = ["severity\tissue_type\tmessage\texpected\tfound"]

    for issue in result.issues:
        expected = str(issue.expected) if issue.expected is not None else ""
        found = str(issue.found) if issue.found is not None else ""
        lines.append(
            f"{issue.severity}\t{issue.issue_type}\t{issue.message}\t{expected}\t{found}"
        )

    return "\n".join(lines)
