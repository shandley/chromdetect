"""
Tests for custom pattern loading functionality.
"""


import pytest

from chromdetect.core import detect_by_name
from chromdetect.patterns import (
    load_custom_patterns,
    merge_patterns,
)


class TestLoadCustomPatterns:
    """Tests for loading custom patterns from files."""

    def test_load_json_patterns(self, tmp_path):
        """Test loading patterns from JSON file."""
        json_content = """{
    "chromosome_patterns": [
        {"pattern": "^MyScaffold_(\\\\d+)$", "name": "my_scaffold"}
    ],
    "unlocalized_patterns": ["my_random"],
    "fragment_patterns": ["my_contig"]
}"""
        json_file = tmp_path / "patterns.json"
        json_file.write_text(json_content)

        chr_patterns, unloc_patterns, frag_patterns = load_custom_patterns(json_file)

        assert len(chr_patterns) == 1
        assert chr_patterns[0][1] == "my_scaffold"
        assert "my_random" in unloc_patterns
        assert "my_contig" in frag_patterns

    def test_load_yaml_patterns(self, tmp_path):
        """Test loading patterns from YAML file."""
        # Try to check if yaml is available - if not, use simpler format
        import importlib.util

        yaml_available = importlib.util.find_spec("yaml") is not None

        if yaml_available:
            # Complex YAML with multiple patterns
            yaml_content = """chromosome_patterns:
  - pattern: "^MyScaffold_(\\\\d+)$"
    name: "my_scaffold"
  - pattern: "^CustomChr_(\\\\d+)$"
    name: "custom_chr"
unlocalized_patterns:
  - my_random
  - my_unloc
fragment_patterns:
  - my_contig
"""
            yaml_file = tmp_path / "patterns.yaml"
            yaml_file.write_text(yaml_content)

            chr_patterns, unloc_patterns, frag_patterns = load_custom_patterns(yaml_file)

            assert len(chr_patterns) == 2
            assert len(unloc_patterns) == 2
            assert len(frag_patterns) == 1
        else:
            # Simple YAML with single pattern (works with basic parser)
            yaml_content = """chromosome_patterns:
  - pattern: "^MyScaffold_(\\\\d+)$"
    name: "my_scaffold"
unlocalized_patterns:
  - my_random
fragment_patterns:
  - my_contig
"""
            yaml_file = tmp_path / "patterns.yaml"
            yaml_file.write_text(yaml_content)

            chr_patterns, unloc_patterns, frag_patterns = load_custom_patterns(yaml_file)

            assert len(chr_patterns) >= 1
            assert len(unloc_patterns) >= 1
            assert len(frag_patterns) >= 1

    def test_load_yaml_minimal(self, tmp_path):
        """Test loading YAML with only chromosome patterns."""
        yaml_content = """chromosome_patterns:
  - pattern: "^MyScaffold_(\\\\d+)$"
    name: "my_scaffold"
"""
        yaml_file = tmp_path / "patterns.yaml"
        yaml_file.write_text(yaml_content)

        chr_patterns, unloc_patterns, frag_patterns = load_custom_patterns(yaml_file)

        assert len(chr_patterns) == 1
        assert len(unloc_patterns) == 0
        assert len(frag_patterns) == 0

    def test_file_not_found(self, tmp_path):
        """Test error for missing file."""
        with pytest.raises(FileNotFoundError):
            load_custom_patterns(tmp_path / "nonexistent.json")

    def test_invalid_json(self, tmp_path):
        """Test error for invalid JSON."""
        json_file = tmp_path / "invalid.json"
        json_file.write_text("{ invalid json }")

        with pytest.raises(ValueError, match="Invalid JSON"):
            load_custom_patterns(json_file)

    def test_invalid_structure(self, tmp_path):
        """Test error for invalid pattern structure."""
        json_content = """{
    "chromosome_patterns": [
        {"pattern": "^Test$"}
    ]
}"""
        json_file = tmp_path / "patterns.json"
        json_file.write_text(json_content)

        with pytest.raises(ValueError, match="pattern.*name"):
            load_custom_patterns(json_file)


class TestMergePatterns:
    """Tests for merging custom patterns with built-in patterns."""

    def test_merge_prepend(self, tmp_path):
        """Test merging patterns with prepend (custom first)."""
        custom_chr = [("^CustomChr_(\\d+)$", "custom_chr")]
        custom_unloc = ["custom_random"]
        custom_frag = ["custom_contig"]

        merged_chr, merged_unloc, merged_frag = merge_patterns(
            custom_chr, custom_unloc, custom_frag, prepend=True
        )

        # Custom patterns should be first
        assert merged_chr[0][1] == "custom_chr"
        # Built-in patterns should still be present
        assert len(merged_chr) > 1

    def test_merge_append(self, tmp_path):
        """Test merging patterns with append (built-in first)."""
        custom_chr = [("^CustomChr_(\\d+)$", "custom_chr")]
        custom_unloc = []
        custom_frag = []

        merged_chr, merged_unloc, merged_frag = merge_patterns(
            custom_chr, custom_unloc, custom_frag, prepend=False
        )

        # Built-in patterns should be first
        assert merged_chr[0][1] != "custom_chr"
        # Custom pattern should be last
        assert merged_chr[-1][1] == "custom_chr"


class TestDetectByNameWithCustomPatterns:
    """Tests for detect_by_name with custom patterns."""

    def test_custom_pattern_detection(self, tmp_path):
        """Test that custom patterns are used for detection."""
        custom_chr = [("^MySpecialScaffold_(\\d+)$", "my_special")]
        custom_unloc = []
        custom_frag = []

        merged = merge_patterns(custom_chr, custom_unloc, custom_frag)

        # Should detect with custom pattern
        result = detect_by_name("MySpecialScaffold_1", custom_patterns=merged)
        classification, confidence, method, chr_id = result

        assert classification == "chromosome"
        assert method == "name_my_special"
        assert chr_id == "1"

    def test_custom_pattern_priority(self, tmp_path):
        """Test that custom patterns have priority when prepended."""
        # Create a custom pattern that overlaps with built-in
        custom_chr = [("^chr(\\d+)$", "custom_chr")]
        custom_unloc = []
        custom_frag = []

        merged = merge_patterns(custom_chr, custom_unloc, custom_frag, prepend=True)

        result = detect_by_name("chr1", custom_patterns=merged)
        classification, confidence, method, chr_id = result

        assert classification == "chromosome"
        # Should use custom pattern since it's first
        assert method == "name_custom_chr"

    def test_custom_unlocalized_pattern(self, tmp_path):
        """Test custom unlocalized pattern detection."""
        custom_chr = []
        custom_unloc = ["special_unloc"]
        custom_frag = []

        merged = merge_patterns(custom_chr, custom_unloc, custom_frag)

        result = detect_by_name("scaffold_special_unloc_001", custom_patterns=merged)
        classification, confidence, method, chr_id = result

        assert classification == "unlocalized"

    def test_custom_fragment_pattern(self, tmp_path):
        """Test custom fragment pattern detection."""
        custom_chr = []
        custom_unloc = []
        custom_frag = ["my_fragment"]

        merged = merge_patterns(custom_chr, custom_unloc, custom_frag)

        result = detect_by_name("scaffold_my_fragment_001", custom_patterns=merged)
        classification, confidence, method, chr_id = result

        assert classification == "unplaced"


class TestPatternFileFormats:
    """Tests for different pattern file formats."""

    def test_yml_extension(self, tmp_path):
        """Test .yml extension is handled."""
        yml_content = """chromosome_patterns:
  - pattern: "^Test_(\\\\d+)$"
    name: "test"
"""
        yml_file = tmp_path / "patterns.yml"
        yml_file.write_text(yml_content)

        chr_patterns, _, _ = load_custom_patterns(yml_file)
        assert len(chr_patterns) == 1

    def test_unknown_extension_json_content(self, tmp_path):
        """Test file with unknown extension but valid JSON content."""
        json_content = '{"chromosome_patterns": [{"pattern": "^Test$", "name": "test"}]}'
        txt_file = tmp_path / "patterns.txt"
        txt_file.write_text(json_content)

        chr_patterns, _, _ = load_custom_patterns(txt_file)
        assert len(chr_patterns) == 1
