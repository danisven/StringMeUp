#!/usr/bin/env python3
"""
Tests for K-mer string processing functions.

CRITICAL: These tests verify the bug fix for ambiguous k-mer handling (v0.2.0).
The bug was that ambiguous k-mers were excluded from the total, inflating confidence scores.
"""

import pytest
import sys
import os

# Add parent directory to path to import stringmeup
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from stringmeup.stringmeup import process_kmer_string


class TestProcessKmerStringBasic:
    """Test basic k-mer string processing without ambiguous k-mers."""

    def test_single_taxon_no_ambiguous(self):
        """Single taxon, no ambiguous k-mers."""
        kmer_string = "562:10"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10}
        assert ambiguous_count == 0

    def test_multiple_taxa_no_ambiguous(self):
        """Multiple taxa, no ambiguous k-mers."""
        kmer_string = "562:10 561:5 0:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 561: 5, 0: 3}
        assert ambiguous_count == 0

    def test_repeated_taxon(self):
        """Same taxon appears multiple times - should aggregate."""
        kmer_string = "562:10 561:5 562:8"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 18, 561: 5}  # 562 appears twice: 10 + 8 = 18
        assert ambiguous_count == 0


class TestProcessKmerStringAmbiguous:
    """CRITICAL TESTS: Verify ambiguous k-mers are handled correctly (v0.2.0 bug fix)."""

    def test_single_ambiguous_segment(self):
        """Single ambiguous segment - MUST be counted separately."""
        kmer_string = "562:10 A:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10}
        assert ambiguous_count == 5  # CRITICAL: Must be 5, not 0

    def test_multiple_ambiguous_segments(self):
        """Multiple ambiguous segments - MUST aggregate all ambiguous counts."""
        kmer_string = "562:10 A:3 561:5 A:7"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 561: 5}
        assert ambiguous_count == 10  # CRITICAL: 3 + 7 = 10

    def test_only_ambiguous(self):
        """Read with only ambiguous k-mers."""
        kmer_string = "A:50"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {}
        assert ambiguous_count == 50

    def test_ambiguous_at_start(self):
        """Ambiguous segment at the start."""
        kmer_string = "A:5 562:10 561:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 561: 3}
        assert ambiguous_count == 5

    def test_ambiguous_at_end(self):
        """Ambiguous segment at the end."""
        kmer_string = "562:10 561:3 A:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 561: 3}
        assert ambiguous_count == 5

    def test_ambiguous_interleaved(self):
        """Ambiguous segments interleaved with taxa."""
        kmer_string = "562:5 A:2 561:3 A:1 0:4 A:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 5, 561: 3, 0: 4}
        assert ambiguous_count == 6  # 2 + 1 + 3 = 6


class TestProcessKmerStringPairedEnd:
    """Test paired-end read processing."""

    def test_paired_end_no_ambiguous(self):
        """Paired-end reads without ambiguous k-mers."""
        kmer_string = "562:10 |:| 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        assert taxa_dict == {562: 10, 561: 5}
        assert ambiguous_count == 0

    def test_paired_end_with_ambiguous_both_reads(self):
        """CRITICAL: Paired-end with ambiguous in BOTH reads."""
        kmer_string = "562:10 A:3 |:| 561:5 A:2"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        assert taxa_dict == {562: 10, 561: 5}
        assert ambiguous_count == 5  # CRITICAL: 3 + 2 = 5 from both reads

    def test_paired_end_with_ambiguous_read1_only(self):
        """Paired-end with ambiguous only in read 1."""
        kmer_string = "562:10 A:5 |:| 561:8"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        assert taxa_dict == {562: 10, 561: 8}
        assert ambiguous_count == 5

    def test_paired_end_with_ambiguous_read2_only(self):
        """Paired-end with ambiguous only in read 2."""
        kmer_string = "562:10 |:| A:5 561:8"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        assert taxa_dict == {562: 10, 561: 8}
        assert ambiguous_count == 5

    def test_paired_end_same_taxon_both_reads(self):
        """Paired-end where same taxon appears in both reads - should aggregate."""
        kmer_string = "562:10 |:| 562:5 561:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        assert taxa_dict == {562: 15, 561: 3}  # 562: 10 + 5 = 15
        assert ambiguous_count == 0

    def test_paired_end_delimiter_not_removed_when_single_end(self):
        """Delimiter should only be removed when paired_input=True."""
        kmer_string = "562:10 |:| 561:5"
        # This is intentionally passing paired_input=False with paired data
        # to test the delimiter handling
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # When not paired mode, |:| is not removed, so it might cause parsing issues
        # The function should handle this gracefully by skipping malformed entries
        assert 562 in taxa_dict or 561 in taxa_dict  # At least some data parsed


class TestProcessKmerStringEdgeCases:
    """Test edge cases and malformed input."""

    def test_empty_string(self):
        """Empty k-mer string."""
        kmer_string = ""
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {}
        assert ambiguous_count == 0

    def test_whitespace_only(self):
        """Only whitespace."""
        kmer_string = "   "
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {}
        assert ambiguous_count == 0

    def test_malformed_no_colon(self):
        """Malformed entry without colon - should skip."""
        kmer_string = "562:10 MALFORMED 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 561: 5}
        assert ambiguous_count == 0

    def test_malformed_non_numeric_count(self):
        """Malformed entry with non-numeric count - should skip."""
        kmer_string = "562:10 561:ABC 0:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 0: 5}
        assert ambiguous_count == 0

    def test_malformed_non_numeric_taxid(self):
        """Malformed entry with non-numeric tax_id - should skip."""
        kmer_string = "562:10 XYZ:5 561:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 561: 3}
        assert ambiguous_count == 0

    def test_reading_frame_borders(self):
        """Reading frame borders (-:-) should be skipped."""
        kmer_string = "562:10 -:- 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 561: 5}
        assert ambiguous_count == 0

    def test_zero_count(self):
        """Entry with zero count."""
        kmer_string = "562:0 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 0, 561: 5}
        assert ambiguous_count == 0

    def test_large_counts(self):
        """Very large k-mer counts."""
        kmer_string = "562:100000 A:50000 561:75000"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 100000, 561: 75000}
        assert ambiguous_count == 50000


class TestProcessKmerStringRegressionBugFix:
    """
    REGRESSION TESTS for v0.2.0 bug fix.

    These tests document the OLD (buggy) behavior vs NEW (correct) behavior.
    The bug was that ambiguous k-mers were excluded from the total.
    """

    def test_bug_fix_ambiguous_counted(self):
        """
        CRITICAL REGRESSION TEST: Verify ambiguous k-mers are now counted.

        Old behavior (v0.x): ambiguous_count would be 0 (filtered out)
        New behavior (v0.2.0): ambiguous_count should be 10
        """
        kmer_string = "562:60 561:20 A:10 0:10"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Old buggy behavior: ambiguous_count == 0
        # New correct behavior: ambiguous_count == 10
        assert ambiguous_count == 10, "BUG FIX VERIFICATION: Ambiguous k-mers MUST be counted"

        # Taxa should not include ambiguous
        assert 'A' not in taxa_dict
        assert taxa_dict == {562: 60, 561: 20, 0: 10}

    def test_bug_fix_return_type(self):
        """
        CRITICAL: Verify function returns tuple, not just dict.

        Old behavior (v0.x): Returned only dict
        New behavior (v0.2.0): Returns (dict, int) tuple
        """
        kmer_string = "562:10 A:5"
        result = process_kmer_string(kmer_string, paired_input=False)

        assert isinstance(result, tuple), "Must return tuple, not just dict"
        assert len(result) == 2, "Must return (taxa_dict, ambiguous_count)"

        taxa_dict, ambiguous_count = result
        assert isinstance(taxa_dict, dict)
        assert isinstance(ambiguous_count, int)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
