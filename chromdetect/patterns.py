"""
Chromosome naming pattern definitions for ChromDetect.

This module contains regex patterns for detecting chromosome-level scaffolds
across different naming conventions used by various genome assembly pipelines.

Adding New Patterns:
    To add support for a new naming convention:
    1. Add the regex pattern to CHROMOSOME_PATTERNS
    2. Include a descriptive method name
    3. Ensure the pattern captures the chromosome ID in group 1
    4. Add tests in tests/test_patterns.py

Custom Patterns:
    You can load custom patterns from a YAML or JSON file using
    load_custom_patterns(). The file format should be:

    YAML:
        chromosome_patterns:
          - pattern: "^MyScaffold_(\\d+)$"
            name: "my_scaffold"
        unlocalized_patterns:
          - "my_random"
        fragment_patterns:
          - "my_contig"

    JSON:
        {
            "chromosome_patterns": [
                {"pattern": "^MyScaffold_(\\d+)$", "name": "my_scaffold"}
            ],
            "unlocalized_patterns": ["my_random"],
            "fragment_patterns": ["my_contig"]
        }
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from re import Pattern

# Common chromosome naming patterns (case-insensitive by default)
# Format: (pattern, method_name)
# The pattern should capture chromosome ID (number, X, Y, etc.) in group 1 when possible
CHROMOSOME_PATTERNS: list[tuple[str, str]] = [
    # Explicit chromosome names (most reliable)
    (r'^chr(?:omosome)?[_\-\s]?(\d+|[XYZWAB]|MT?|Un)$', 'chr_explicit'),

    # VGP-style naming conventions
    (r'^super[_\-\s]?scaffold[_\-\s]?(\d+|[XYZWAB])$', 'super_scaffold'),
    (r'^superscaffold[_\-\s]?(\d+|[XYZWAB])$', 'superscaffold'),
    (r'^SUPER[_\-\s]?(\d+|[XYZWAB])$', 'SUPER'),

    # Linkage group naming (common in fish, plants)
    (r'^LG[_\-\s]?(\d+|[XYZWAB])$', 'linkage_group'),

    # GenBank/RefSeq accession patterns (assembled chromosomes)
    (r'^NC_\d+\.\d+$', 'ncbi_refseq'),
    (r'^CM\d+\.\d+$', 'ncbi_genbank'),

    # Simple numeric chromosome names
    (r'^(?:chr)?(\d+|[XYZWAB])$', 'numeric'),

    # HiC scaffolder output patterns
    (r'^HiC_scaffold_(\d+)$', 'hic_scaffold'),
    (r'^Scaffold_(\d+)_RaGOO$', 'ragoo'),

    # Plant genome conventions
    (r'^Pt(\d+)$', 'plant_chromosome'),
    (r'^Gm(\d+)$', 'soybean_chromosome'),

    # Additional VGP variations
    (r'^Super-Scaffold_(\d+)$', 'super_scaffold_hyphen'),
    (r'^scaffold_(\d+)_cov\d+$', 'scaffold_cov'),
]

# Patterns suggesting unlocalized/random scaffolds (chromosome-associated but not placed)
UNLOCALIZED_PATTERNS: list[str] = [
    r'random',
    r'unloc',
    r'unplaced',
    r'_un_',
    r'chrUn',
    r'scaffold.*unloc',
    r'_hap\d+_unloc',
]

# Patterns suggesting contigs/fragments (not chromosome-level)
FRAGMENT_PATTERNS: list[str] = [
    r'ctg\d*$',
    r'contig',
    r'_arrow_',
    r'_pilon',
    r'fragment',
    r'_hap\d',
    r'_alt$',
    r'_patch$',
    r'debris',
]

# Sex chromosome identifiers
SEX_CHROMOSOMES: set[str] = {'X', 'Y', 'Z', 'W', 'x', 'y', 'z', 'w'}

# Mitochondrial identifiers
MITOCHONDRIAL_IDS: set[str] = {'M', 'MT', 'Mt', 'mt', 'mito', 'mitochondrion'}


def compile_patterns() -> list[tuple[Pattern[str], str]]:
    """Compile chromosome patterns for efficient matching."""
    return [(re.compile(p, re.IGNORECASE), name) for p, name in CHROMOSOME_PATTERNS]


def compile_exclusion_patterns() -> tuple[list[Pattern[str]], list[Pattern[str]]]:
    """Compile exclusion patterns (unlocalized, fragments)."""
    unloc = [re.compile(p, re.IGNORECASE) for p in UNLOCALIZED_PATTERNS]
    frag = [re.compile(p, re.IGNORECASE) for p in FRAGMENT_PATTERNS]
    return unloc, frag


# Pre-compiled patterns for performance
COMPILED_CHROMOSOME_PATTERNS = compile_patterns()
COMPILED_UNLOCALIZED, COMPILED_FRAGMENT = compile_exclusion_patterns()


def load_custom_patterns(
    file_path: Path | str,
) -> tuple[
    list[tuple[str, str]],
    list[str],
    list[str],
]:
    """
    Load custom patterns from a YAML or JSON file.

    The file can be either YAML or JSON format. The structure should be:
    - chromosome_patterns: list of {pattern: str, name: str} objects
    - unlocalized_patterns: list of pattern strings (optional)
    - fragment_patterns: list of pattern strings (optional)

    Args:
        file_path: Path to the patterns file (YAML or JSON)

    Returns:
        Tuple of (chromosome_patterns, unlocalized_patterns, fragment_patterns)

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Patterns file not found: {file_path}")

    # Read file content
    content = file_path.read_text()

    # Try to parse as JSON first, then YAML
    data = None
    if file_path.suffix.lower() in (".json",):
        try:
            data = json.loads(content)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON format: {e}") from e
    elif file_path.suffix.lower() in (".yaml", ".yml"):
        # Try to import yaml, fall back to basic parsing
        try:
            import yaml  # type: ignore[import-untyped]

            data = yaml.safe_load(content)
        except ImportError:
            # Basic YAML parsing for simple structures
            data = _parse_simple_yaml(content)
        except Exception as e:
            raise ValueError(f"Invalid YAML format: {e}") from e
    else:
        # Try JSON first, then YAML
        try:
            data = json.loads(content)
        except json.JSONDecodeError:
            try:
                import yaml  # type: ignore[import-untyped]

                data = yaml.safe_load(content)
            except ImportError:
                data = _parse_simple_yaml(content)
            except Exception as e:
                raise ValueError(f"Could not parse file as JSON or YAML: {e}") from e

    if not isinstance(data, dict):
        raise ValueError("Patterns file must contain a dictionary/object at root level")

    # Extract patterns
    chr_patterns: list[tuple[str, str]] = []
    if "chromosome_patterns" in data:
        for item in data["chromosome_patterns"]:
            if isinstance(item, dict) and "pattern" in item and "name" in item:
                chr_patterns.append((item["pattern"], item["name"]))
            else:
                raise ValueError(
                    "chromosome_patterns must be list of {pattern, name} objects"
                )

    unloc_patterns: list[str] = data.get("unlocalized_patterns", [])
    frag_patterns: list[str] = data.get("fragment_patterns", [])

    return chr_patterns, unloc_patterns, frag_patterns


def _parse_simple_yaml(content: str) -> dict:
    """
    Very basic YAML parser for simple key-value and list structures.

    This is a fallback when PyYAML is not installed.
    Only handles simple cases like:
        key:
          - item1
          - item2
        nested:
          - pattern: "..."
            name: "..."
    """
    lines = content.strip().split("\n")
    result: dict = {}
    current_key: str | None = None
    current_list: list = []
    current_dict: dict = {}

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        # Calculate indent
        line_indent = len(line) - len(line.lstrip())

        # Top-level key
        if line_indent == 0 and ":" in stripped:
            # Save previous list if any
            if current_key and current_list:
                result[current_key] = current_list
                current_list = []
            if current_key and current_dict:
                current_list.append(current_dict)
                result[current_key] = current_list
                current_list = []
                current_dict = {}

            key = stripped.split(":")[0].strip()
            value = ":".join(stripped.split(":")[1:]).strip()
            current_key = key
            if value:
                result[key] = value

        # List item
        elif stripped.startswith("-"):
            item_content = stripped[1:].strip()
            # Check if it's a key-value in list item
            if ":" in item_content and not item_content.startswith('"'):
                if current_dict:
                    current_list.append(current_dict)
                current_dict = {}
                key, value = item_content.split(":", 1)
                current_dict[key.strip()] = value.strip().strip('"').strip("'")
            else:
                if current_dict:
                    current_list.append(current_dict)
                    current_dict = {}
                current_list.append(item_content.strip('"').strip("'"))

        # Continuation of dict item
        elif ":" in stripped and current_dict is not None:
            key, value = stripped.split(":", 1)
            current_dict[key.strip()] = value.strip().strip('"').strip("'")

    # Save final items
    if current_dict:
        current_list.append(current_dict)
    if current_key and current_list:
        result[current_key] = current_list

    return result


def merge_patterns(
    custom_chr: list[tuple[str, str]],
    custom_unloc: list[str],
    custom_frag: list[str],
    prepend: bool = True,
) -> tuple[
    list[tuple[Pattern[str], str]],
    list[Pattern[str]],
    list[Pattern[str]],
]:
    """
    Merge custom patterns with built-in patterns.

    Args:
        custom_chr: Custom chromosome patterns
        custom_unloc: Custom unlocalized patterns
        custom_frag: Custom fragment patterns
        prepend: If True, custom patterns are checked first (higher priority)

    Returns:
        Tuple of compiled (chromosome_patterns, unlocalized_patterns, fragment_patterns)
    """
    # Compile custom patterns
    compiled_custom_chr = [(re.compile(p, re.IGNORECASE), n) for p, n in custom_chr]
    compiled_custom_unloc = [re.compile(p, re.IGNORECASE) for p in custom_unloc]
    compiled_custom_frag = [re.compile(p, re.IGNORECASE) for p in custom_frag]

    if prepend:
        merged_chr = compiled_custom_chr + COMPILED_CHROMOSOME_PATTERNS
        merged_unloc = compiled_custom_unloc + COMPILED_UNLOCALIZED
        merged_frag = compiled_custom_frag + COMPILED_FRAGMENT
    else:
        merged_chr = COMPILED_CHROMOSOME_PATTERNS + compiled_custom_chr
        merged_unloc = COMPILED_UNLOCALIZED + compiled_custom_unloc
        merged_frag = COMPILED_FRAGMENT + compiled_custom_frag

    return merged_chr, merged_unloc, merged_frag
