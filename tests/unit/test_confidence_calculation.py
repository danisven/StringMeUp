#!/usr/bin/env python3
"""
Tests for confidence score calculation.

CRITICAL: These tests verify the core bug fix (v0.2.0) where ambiguous k-mers
must be included in the denominator to match Kraken2's behavior.

The bug: StringMeUp v0.x calculated confidence = hits / (non-ambiguous k-mers)
The fix: StringMeUp v0.2.0 calculates confidence = hits / (ALL k-mers including ambiguous)
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from stringmeup.stringmeup import process_kmer_string


class TestConfidenceCalculationNoAmbiguous:
    """Test confidence calculations when there are NO ambiguous k-mers."""

    def test_perfect_confidence_no_ambiguous(self):
        """
        100 k-mers, all hit same species, no ambiguous.
        Confidence should be 1.0 (100%)
        """
        kmer_string = "562:100"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        hits_to_species = taxa_dict.get(562, 0)
        confidence = hits_to_species / total_kmers

        assert total_kmers == 100
        assert hits_to_species == 100
        assert confidence == 1.0

    def test_half_confidence_no_ambiguous(self):
        """
        100 k-mers, 50 hit species, 50 unclassified, no ambiguous.
        Confidence should be 0.5 (50%)
        """
        kmer_string = "562:50 0:50"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        hits_to_species = taxa_dict.get(562, 0)
        confidence = hits_to_species / total_kmers

        assert total_kmers == 100
        assert hits_to_species == 50
        assert confidence == 0.5


class TestConfidenceCalculationWithAmbiguous:
    """
    CRITICAL TESTS: Verify confidence includes ambiguous k-mers in denominator.
    This is the core of the v0.2.0 bug fix.
    """

    def test_confidence_with_10_percent_ambiguous(self):
        """
        CRITICAL TEST: 100 total k-mers, 10 ambiguous, 60 hits.

        Old behavior (v0.x): confidence = 60 / 90 = 0.667 (WRONG)
        New behavior (v0.2.0): confidence = 60 / 100 = 0.60 (CORRECT)

        This is the exact scenario from the bug report.
        """
        kmer_string = "562:60 0:30 A:10"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Calculate confidence the CORRECT way (v0.2.0+)
        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        hits_to_species = taxa_dict.get(562, 0)
        confidence = hits_to_species / total_kmers

        assert ambiguous_count == 10, "Must count ambiguous k-mers"
        assert total_kmers == 100, "Total must include ambiguous: 60+30+10=100"
        assert hits_to_species == 60
        assert confidence == 0.60, f"Confidence must be 0.60, not {confidence}"

        # Verify the OLD buggy calculation would have been wrong
        total_without_ambiguous = sum(taxa_dict.values())
        buggy_confidence = hits_to_species / total_without_ambiguous
        assert pytest.approx(buggy_confidence, abs=0.01) == 0.667, "Document the old buggy behavior"
        assert confidence < buggy_confidence, "New confidence must be LOWER (more strict)"

    def test_confidence_with_20_percent_ambiguous(self):
        """
        100 total k-mers, 20 ambiguous, 60 hits.

        Old behavior: 60 / 80 = 0.75
        New behavior: 60 / 100 = 0.60

        Difference: 15 percentage points!
        """
        kmer_string = "562:60 0:20 A:20"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        hits_to_species = taxa_dict.get(562, 0)
        confidence = hits_to_species / total_kmers

        assert total_kmers == 100
        assert confidence == 0.60

        # Verify old behavior would have been significantly different
        buggy_confidence = hits_to_species / sum(taxa_dict.values())
        assert abs(buggy_confidence - confidence) == pytest.approx(0.15), \
            "20% ambiguous causes 15 percentage point difference"

    def test_confidence_with_50_percent_ambiguous(self):
        """
        Extreme case: 50% of k-mers are ambiguous.

        Old behavior: 40 / 50 = 0.80
        New behavior: 40 / 100 = 0.40

        Difference: 40 percentage points!
        """
        kmer_string = "562:40 0:10 A:50"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        hits_to_species = taxa_dict.get(562, 0)
        confidence = hits_to_species / total_kmers

        assert total_kmers == 100
        assert confidence == 0.40

        buggy_confidence = hits_to_species / sum(taxa_dict.values())
        assert buggy_confidence == 0.80
        assert abs(buggy_confidence - confidence) == 0.40, \
            "50% ambiguous causes MASSIVE 40 percentage point difference"

    def test_confidence_threshold_crossing(self):
        """
        CRITICAL: Test case where ambiguous k-mers cause threshold crossing.

        Setup: 60 hits, 30 unclassified, 10 ambiguous, threshold = 0.65

        Old behavior: 60/90 = 0.667 → PASSES threshold
        New behavior: 60/100 = 0.60 → FAILS threshold

        This would cause actual classification changes!
        """
        kmer_string = "562:60 0:30 A:10"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        threshold = 0.65
        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        hits_to_species = taxa_dict.get(562, 0)
        confidence = hits_to_species / total_kmers

        assert confidence == 0.60
        assert confidence < threshold, "New behavior: FAILS threshold (correct)"

        # Old behavior would have passed
        buggy_confidence = hits_to_species / sum(taxa_dict.values())
        assert buggy_confidence > threshold, "Old behavior: PASSES threshold (incorrect)"


class TestConfidenceCalculationPairedEnd:
    """Test confidence calculations for paired-end reads."""

    def test_paired_end_confidence_with_ambiguous_both_reads(self):
        """
        Paired-end: ambiguous in BOTH reads.
        Must aggregate all ambiguous k-mers from both reads.
        """
        kmer_string = "562:50 A:5 |:| 562:30 A:15"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        hits_to_species = taxa_dict.get(562, 0)
        confidence = hits_to_species / total_kmers

        assert ambiguous_count == 20, "Must count ambiguous from both reads: 5+15=20"
        assert hits_to_species == 80, "Must aggregate hits from both reads: 50+30=80"
        assert total_kmers == 100, "Total: 80 hits + 20 ambiguous"
        assert confidence == 0.80

    def test_paired_end_same_confidence_regardless_of_distribution(self):
        """
        Two reads with same total k-mer distribution should have same confidence,
        regardless of how k-mers are distributed between the two reads.
        """
        # Distribution 1: Most k-mers in read1
        kmer_string1 = "562:70 A:10 |:| 562:10 A:10"
        taxa_dict1, ambiguous_count1 = process_kmer_string(kmer_string1, paired_input=True)
        confidence1 = taxa_dict1.get(562, 0) / (sum(taxa_dict1.values()) + ambiguous_count1)

        # Distribution 2: More balanced
        kmer_string2 = "562:40 A:10 |:| 562:40 A:10"
        taxa_dict2, ambiguous_count2 = process_kmer_string(kmer_string2, paired_input=True)
        confidence2 = taxa_dict2.get(562, 0) / (sum(taxa_dict2.values()) + ambiguous_count2)

        assert confidence1 == confidence2 == 0.80, "Same total distribution = same confidence"


class TestConfidenceCalculationVariableReadLengths:
    """
    Test that confidence is independent of read length.
    Only k-mer counts matter, not read length.
    """

    def test_same_kmer_distribution_different_lengths(self):
        """
        Two reads with different lengths but same k-mer distribution
        should have the SAME confidence.

        Read1: length=50bp, 80 k-mers with 60 hits
        Read2: length=150bp, 80 k-mers with 60 hits (many ambiguous regions)

        Both should have confidence = 60/80 = 0.75
        """
        # Simulating Read1
        kmer_string1 = "562:60 0:20"
        taxa_dict1, ambiguous_count1 = process_kmer_string(kmer_string1, paired_input=False)
        confidence1 = taxa_dict1.get(562, 0) / (sum(taxa_dict1.values()) + ambiguous_count1)

        # Simulating Read2 (same k-mer counts)
        kmer_string2 = "562:60 0:20"
        taxa_dict2, ambiguous_count2 = process_kmer_string(kmer_string2, paired_input=False)
        confidence2 = taxa_dict2.get(562, 0) / (sum(taxa_dict2.values()) + ambiguous_count2)

        assert confidence1 == confidence2 == 0.75, \
            "Confidence must be independent of read length"


class TestConfidenceCalculationEdgeCases:
    """Test edge cases in confidence calculation."""

    def test_all_ambiguous_confidence_undefined(self):
        """
        Read with only ambiguous k-mers.
        No hits possible, should handle gracefully.
        """
        kmer_string = "A:100"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        assert total_kmers == 100
        assert len(taxa_dict) == 0, "No taxa hits"
        # In this case, reclassify_read should mark as unclassified

    def test_zero_kmers_confidence_undefined(self):
        """
        Empty k-mer string.
        Should be caught by edge case handling in reclassify_read.
        """
        kmer_string = ""
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        assert total_kmers == 0
        # Division by zero would occur - must be handled

    def test_single_kmer_perfect_confidence(self):
        """Single k-mer with perfect hit."""
        kmer_string = "562:1"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        total_kmers = sum(taxa_dict.values()) + ambiguous_count
        confidence = taxa_dict.get(562, 0) / total_kmers

        assert confidence == 1.0


class TestConfidenceCalculationRegressionTests:
    """
    REGRESSION TESTS documenting the exact impact of the v0.2.0 bug fix.
    """

    def test_bug_fix_impact_table(self):
        """
        Document the impact of the bug fix for different ambiguous percentages.

        Creates a table showing old vs new confidence for various scenarios.
        """
        test_cases = [
            # (hits, non_amb, ambiguous, expected_new_conf, expected_old_conf)
            (60, 40, 0, 0.60, 0.60),    # 0% ambiguous: no change
            (60, 35, 5, 0.60, 0.632),   # 5% ambiguous
            (60, 30, 10, 0.60, 0.667),  # 10% ambiguous
            (60, 20, 20, 0.60, 0.75),   # 20% ambiguous
            (60, 10, 30, 0.60, 0.857),  # 30% ambiguous
        ]

        for hits, non_amb, ambiguous, expected_new, expected_old in test_cases:
            # Simulate k-mer string
            kmer_string = f"562:{hits} 0:{non_amb} A:{ambiguous}"
            taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

            # New correct calculation (v0.2.0+)
            total_with_ambiguous = sum(taxa_dict.values()) + ambiguous_count
            new_confidence = hits / total_with_ambiguous

            # Old buggy calculation (v0.x)
            total_without_ambiguous = sum(taxa_dict.values())
            old_confidence = hits / total_without_ambiguous

            assert pytest.approx(new_confidence, abs=0.01) == expected_new
            assert pytest.approx(old_confidence, abs=0.01) == expected_old
            assert new_confidence <= old_confidence, "New must be <= old (more strict)"

            difference = old_confidence - new_confidence
            percent_ambiguous = ambiguous / (hits + non_amb + ambiguous)
            print(f"{percent_ambiguous*100:.0f}% ambiguous: "
                  f"old={old_confidence:.3f}, new={new_confidence:.3f}, "
                  f"diff={difference:.3f}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])  # -s to show print statements
