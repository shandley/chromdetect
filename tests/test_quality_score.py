"""
Tests for quality score calculation functionality.
"""


from chromdetect.core import (
    ScaffoldInfo,
    calculate_quality_score,
)


class TestCalculateQualityScore:
    """Tests for the calculate_quality_score function."""

    def test_quality_score_empty(self):
        """Test quality score with empty results."""
        score = calculate_quality_score([])
        assert score == 0.0

    def test_quality_score_no_chromosomes(self):
        """Test quality score with no chromosomes."""
        results = [
            ScaffoldInfo(
                name="scaffold1",
                length=1000000,
                classification="unplaced",
                confidence=0.5,
                detection_method="name_none",
            )
        ]
        score = calculate_quality_score(results)
        assert score == 0.0

    def test_quality_score_high_confidence(self):
        """Test quality score with high confidence chromosomes."""
        results = [
            ScaffoldInfo(
                name="chr1",
                length=100000000,
                classification="chromosome",
                confidence=0.95,
                detection_method="name_chr_explicit",
                chromosome_id="1",
            ),
            ScaffoldInfo(
                name="chr2",
                length=90000000,
                classification="chromosome",
                confidence=0.95,
                detection_method="name_chr_explicit",
                chromosome_id="2",
            ),
        ]
        score = calculate_quality_score(results)
        # High confidence should contribute positively
        assert score > 0.3

    def test_quality_score_with_telomeres(self):
        """Test quality score with telomere detection."""
        # Without telomeres
        results_no_telo = [
            ScaffoldInfo(
                name="chr1",
                length=100000000,
                classification="chromosome",
                confidence=0.9,
                detection_method="name_chr_explicit",
                has_telomere_5prime=False,
                has_telomere_3prime=False,
            ),
        ]
        score_no_telo = calculate_quality_score(results_no_telo)

        # With T2T chromosomes
        results_t2t = [
            ScaffoldInfo(
                name="chr1",
                length=100000000,
                classification="chromosome",
                confidence=0.9,
                detection_method="name_chr_explicit",
                has_telomere_5prime=True,
                has_telomere_3prime=True,
            ),
        ]
        score_t2t = calculate_quality_score(results_t2t)

        # T2T should have higher score
        assert score_t2t > score_no_telo

    def test_quality_score_with_expected_chromosomes(self):
        """Test quality score with expected chromosome count."""
        results = [
            ScaffoldInfo(
                name="chr1",
                length=100000000,
                classification="chromosome",
                confidence=0.9,
                detection_method="name_chr_explicit",
            ),
            ScaffoldInfo(
                name="chr2",
                length=90000000,
                classification="chromosome",
                confidence=0.9,
                detection_method="name_chr_explicit",
            ),
        ]

        # Exact match
        score_exact = calculate_quality_score(results, expected_chromosomes=2)

        # Too few chromosomes found
        score_fewer = calculate_quality_score(results, expected_chromosomes=10)

        # Exact match should score higher on completeness
        assert score_exact > score_fewer

    def test_quality_score_overclassification(self):
        """Test quality score penalty for over-classification."""
        results = [
            ScaffoldInfo(
                name=f"chr{i}",
                length=10000000,
                classification="chromosome",
                confidence=0.9,
                detection_method="name_chr_explicit",
            )
            for i in range(50)
        ]

        # Expected 10, found 50 - should be penalized
        score = calculate_quality_score(results, expected_chromosomes=10)

        # Should still be positive but lower than if we had 10
        assert 0 < score < 0.8

    def test_quality_score_size_consistency(self):
        """Test quality score with size consistency factor."""
        chr_length = 100000000
        total_length = 110000000  # Chromosomes are 90% of total

        results = [
            ScaffoldInfo(
                name="chr1",
                length=chr_length,
                classification="chromosome",
                confidence=0.9,
                detection_method="name_chr_explicit",
            ),
        ]

        score = calculate_quality_score(results, total_length=total_length)
        # Good size ratio should contribute positively
        assert score > 0.4

    def test_quality_score_bounds(self):
        """Test that quality score is always between 0 and 1."""
        # Very high quality assembly
        results = [
            ScaffoldInfo(
                name=f"chr{i}",
                length=100000000,
                classification="chromosome",
                confidence=0.99,
                detection_method="name_chr_explicit",
                has_telomere_5prime=True,
                has_telomere_3prime=True,
            )
            for i in range(23)
        ]

        score = calculate_quality_score(
            results,
            expected_chromosomes=23,
            total_length=2300000000,
        )

        assert 0.0 <= score <= 1.0

    def test_quality_score_reasonable_count_without_expected(self):
        """Test quality score with reasonable chromosome count but no expected."""
        # 23 chromosomes without expected count
        results = [
            ScaffoldInfo(
                name=f"chr{i}",
                length=100000000,
                classification="chromosome",
                confidence=0.9,
                detection_method="name_chr_explicit",
            )
            for i in range(23)
        ]

        score = calculate_quality_score(results)
        # Should get reasonable completeness score for 10-100 range
        assert score > 0.3
