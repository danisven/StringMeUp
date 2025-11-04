#!/usr/bin/env python3
"""
Tests for helper functions in stringmeup.py.

Tests the utility functions that support the main reclassification logic.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from stringmeup.stringmeup import (
    _calculate_max_confidence,
    _check_early_termination,
)


class TestCalculateMaxConfidence:
    """Tests for _calculate_max_confidence() helper function."""

    def test_max_confidence_single_taxon_no_ambiguous(self):
        """
        Single taxon, no ambiguous k-mers.
        Max confidence should be 1.0 (100%).
        """
        taxa_dict = {562: 100}
        total_kmers = 100

        max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        assert max_conf == 1.0

    def test_max_confidence_single_taxon_with_ambiguous(self):
        """
        CRITICAL: Single taxon with ambiguous k-mers.
        60 hits, 30 unclassified, 10 ambiguous = 100 total.
        Max confidence = 60/100 = 0.60 (NOT 60/90 = 0.667).
        """
        taxa_dict = {562: 60, 0: 30}
        total_kmers = 100  # Includes 10 ambiguous

        max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        assert max_conf == 0.60, "Must include ambiguous in denominator"

    def test_max_confidence_excludes_unclassified(self):
        """
        Max confidence excludes tax_id 0 (unclassified) from numerator.
        """
        taxa_dict = {562: 50, 0: 40}
        total_kmers = 100  # 50 + 40 + 10 ambiguous

        max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        # Should be 50/100, not (50+40)/100
        assert max_conf == 0.50

    def test_max_confidence_multiple_taxa(self):
        """
        Multiple taxa: sum all hits excluding unclassified.
        """
        taxa_dict = {562: 30, 590: 20, 2188: 10, 0: 40}
        total_kmers = 100

        max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        # 30 + 20 + 10 = 60 hits, 100 total
        assert max_conf == 0.60

    def test_max_confidence_only_unclassified(self):
        """
        Only unclassified k-mers (0) and ambiguous.
        Max confidence should be 0.0.
        """
        taxa_dict = {0: 90}
        total_kmers = 100  # 90 + 10 ambiguous

        max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        assert max_conf == 0.0

    def test_max_confidence_only_ambiguous(self):
        """
        Only ambiguous k-mers, no hits.
        Max confidence should be 0.0.
        """
        taxa_dict = {}
        total_kmers = 100  # All ambiguous

        max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        assert max_conf == 0.0

    def test_max_confidence_very_small_values(self):
        """Test with single k-mer."""
        taxa_dict = {562: 1}
        total_kmers = 1

        max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        assert max_conf == 1.0

    def test_max_confidence_large_values(self):
        """Test with large k-mer counts."""
        taxa_dict = {562: 1000000, 590: 500000, 0: 500000}
        total_kmers = 2000000

        max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        # (1000000 + 500000) / 2000000 = 0.75
        assert max_conf == 0.75


class TestCheckEarlyTermination:
    """Tests for _check_early_termination() helper function."""

    def test_should_not_terminate_max_conf_above_threshold(self):
        """
        Max confidence above threshold, sufficient hit groups.
        Should NOT terminate early.
        """
        should_terminate = _check_early_termination(
            max_confidence=0.80,
            confidence_threshold=0.65,
            verbose_input=False,
            minimizer_hit_groups=5,
            minimum_hit_groups=2
        )

        assert should_terminate == False

    def test_should_terminate_max_conf_below_threshold(self):
        """
        Max confidence below threshold.
        Should terminate early (doomed to fail).
        """
        should_terminate = _check_early_termination(
            max_confidence=0.50,
            confidence_threshold=0.65,
            verbose_input=False,
            minimizer_hit_groups=5,
            minimum_hit_groups=2
        )

        assert should_terminate == True

    def test_should_terminate_insufficient_hit_groups(self):
        """
        Sufficient confidence but insufficient hit groups (with verbose_input=True).
        Should terminate early.
        """
        should_terminate = _check_early_termination(
            max_confidence=0.80,
            confidence_threshold=0.65,
            verbose_input=True,  # Must be True for hit groups check
            minimizer_hit_groups=1,  # Below minimum
            minimum_hit_groups=2
        )

        assert should_terminate == True

    def test_threshold_crossing_exactly_at_threshold(self):
        """
        Max confidence exactly at threshold.
        Should NOT terminate (>= threshold passes).
        """
        should_terminate = _check_early_termination(
            max_confidence=0.65,
            confidence_threshold=0.65,
            verbose_input=False,
            minimizer_hit_groups=5,
            minimum_hit_groups=2
        )

        assert should_terminate == False

    def test_hit_groups_exactly_at_minimum(self):
        """
        Hit groups exactly at minimum (with verbose_input=True).
        Should NOT terminate (>= minimum passes).
        """
        should_terminate = _check_early_termination(
            max_confidence=0.80,
            confidence_threshold=0.65,
            verbose_input=True,  # Must be True for hit groups check
            minimizer_hit_groups=2,  # Exactly at minimum
            minimum_hit_groups=2
        )

        assert should_terminate == False

    def test_hit_groups_ignored_when_not_verbose(self):
        """
        Hit groups below minimum but verbose_input=False.
        Should NOT terminate (hit groups only checked when verbose).
        """
        should_terminate = _check_early_termination(
            max_confidence=0.80,
            confidence_threshold=0.65,
            verbose_input=False,  # Hit groups check is skipped
            minimizer_hit_groups=1,  # Below minimum
            minimum_hit_groups=2
        )

        assert should_terminate == False

    def test_both_conditions_fail(self):
        """
        Both conditions fail (low confidence AND low hit groups with verbose=True).
        Should terminate early (either condition sufficient).
        """
        should_terminate = _check_early_termination(
            max_confidence=0.50,
            confidence_threshold=0.65,
            verbose_input=True,  # Enable hit groups check
            minimizer_hit_groups=1,
            minimum_hit_groups=2
        )

        assert should_terminate == True

    def test_zero_confidence(self):
        """
        Zero max confidence.
        Should terminate early.
        """
        should_terminate = _check_early_termination(
            max_confidence=0.0,
            confidence_threshold=0.65,
            verbose_input=False,
            minimizer_hit_groups=5,
            minimum_hit_groups=2
        )

        assert should_terminate == True

    def test_zero_threshold(self):
        """
        Zero threshold (accept anything).
        Should NOT terminate if hit groups sufficient.
        """
        should_terminate = _check_early_termination(
            max_confidence=0.50,
            confidence_threshold=0.0,
            verbose_input=False,
            minimizer_hit_groups=5,
            minimum_hit_groups=2
        )

        assert should_terminate == False


class TestConfidenceFormulaRegressionTests:
    """
    Regression tests documenting the impact of including ambiguous k-mers.
    These tests verify that max_confidence correctly uses total k-mers
    (including ambiguous) in the denominator.
    """

    def test_10_percent_ambiguous_impact(self):
        """
        10% ambiguous k-mers: old vs new max_confidence.

        Old (buggy): 60 / 90 = 0.667
        New (correct): 60 / 100 = 0.60
        """
        taxa_dict = {562: 60, 0: 30}
        total_kmers = 100  # Includes 10 ambiguous

        new_max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        # Calculate what the old buggy version would have given
        old_total = sum(taxa_dict.values())  # 90 (excluded ambiguous)
        old_max_conf = 60 / old_total

        assert new_max_conf == 0.60
        assert pytest.approx(old_max_conf, abs=0.01) == 0.667
        assert new_max_conf < old_max_conf, "New formula more strict"

    def test_20_percent_ambiguous_impact(self):
        """
        20% ambiguous k-mers: shows larger impact.

        Old (buggy): 60 / 80 = 0.75
        New (correct): 60 / 100 = 0.60
        """
        taxa_dict = {562: 60, 0: 20}
        total_kmers = 100  # Includes 20 ambiguous

        new_max_conf = _calculate_max_confidence(taxa_dict, total_kmers)

        old_total = sum(taxa_dict.values())  # 80
        old_max_conf = 60 / old_total

        assert new_max_conf == 0.60
        assert old_max_conf == 0.75
        assert pytest.approx(abs(old_max_conf - new_max_conf), abs=0.001) == 0.15, "15 percentage point difference"

    def test_threshold_crossing_due_to_ambiguous(self):
        """
        CRITICAL: Test where ambiguous k-mers cause threshold crossing.

        Setup: 60 hits, 30 unclassified, 10 ambiguous, threshold = 0.65

        Old: 60/90 = 0.667 → Would NOT terminate (passes threshold)
        New: 60/100 = 0.60 → WILL terminate (fails threshold)
        """
        taxa_dict = {562: 60, 0: 30}
        total_kmers = 100
        threshold = 0.65

        new_max_conf = _calculate_max_confidence(taxa_dict, total_kmers)
        new_terminates = _check_early_termination(
            new_max_conf, threshold, False, 5, 2
        )

        # Old behavior
        old_total = sum(taxa_dict.values())
        old_max_conf = 60 / old_total
        old_terminates = _check_early_termination(
            old_max_conf, threshold, False, 5, 2
        )

        assert new_max_conf == 0.60
        assert new_terminates == True, "New: should terminate (fails threshold)"
        assert old_max_conf > threshold
        assert old_terminates == False, "Old: would not terminate (passes threshold)"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
