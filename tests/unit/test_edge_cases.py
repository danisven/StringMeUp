#!/usr/bin/env python3
"""
Edge case and boundary condition tests.

Tests unusual inputs, error conditions, and boundary values.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from stringmeup.stringmeup import process_kmer_string, ReadClassification


class TestBoundaryValues:
    """Test boundary values and extreme cases."""

    def test_single_kmer(self):
        """Minimum case: single k-mer."""
        kmer_string = "562:1"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 1}
        assert ambiguous_count == 0

    def test_very_large_kmer_count(self):
        """Maximum case: very large k-mer count."""
        kmer_string = "562:1000000"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 1000000}
        assert ambiguous_count == 0

    def test_many_different_taxa(self):
        """Many different taxa in single read."""
        # Create string with 100 different taxa
        taxa_strings = [f"{i}:1" for i in range(100)]
        kmer_string = " ".join(taxa_strings)

        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert len(taxa_dict) == 100
        assert sum(taxa_dict.values()) == 100
        assert ambiguous_count == 0

    def test_very_long_kmer_string(self):
        """Very long k-mer string with many segments."""
        # 1000 segments alternating between taxa and ambiguous
        segments = []
        for i in range(500):
            segments.append(f"562:{i+1}")
            segments.append(f"A:{i+1}")

        kmer_string = " ".join(segments)
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        expected_taxa_total = sum(range(1, 501))  # 1+2+...+500
        expected_ambiguous_total = sum(range(1, 501))

        assert taxa_dict[562] == expected_taxa_total
        assert ambiguous_count == expected_ambiguous_total


class TestMalformedInput:
    """Test handling of malformed and invalid input."""

    def test_completely_invalid_format(self):
        """Completely invalid k-mer string."""
        kmer_string = "this is not a valid kmer string at all"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Should return empty dict and 0 ambiguous (all entries skipped)
        assert taxa_dict == {}
        assert ambiguous_count == 0

    def test_mix_valid_and_invalid(self):
        """Mix of valid and invalid entries."""
        kmer_string = "562:10 INVALID ABC:XYZ 561:5 !@#$ A:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Should process only valid entries
        assert taxa_dict == {562: 10, 561: 5}
        assert ambiguous_count == 3

    def test_negative_counts(self):
        """Negative k-mer counts (should be skipped or handled)."""
        kmer_string = "562:-10 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Negative count parsed as-is (Python int() accepts negative)
        # This is edge case behavior - might want to add validation
        assert 562 in taxa_dict
        assert taxa_dict[561] == 5

    def test_float_instead_of_int(self):
        """Float values instead of integers."""
        kmer_string = "562:10.5 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # int("10.5") raises ValueError, so entry should be skipped
        assert 562 not in taxa_dict
        assert taxa_dict == {561: 5}

    def test_missing_colon(self):
        """Entries missing the colon separator."""
        kmer_string = "56210 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Entry without colon should be skipped
        assert taxa_dict == {561: 5}

    def test_extra_colons(self):
        """Entries with extra colons."""
        kmer_string = "562:10:20 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Will split into 3 parts, len(parts) != 2, should skip
        assert taxa_dict == {561: 5}

    def test_only_colon(self):
        """Entry that is only a colon."""
        kmer_string = ": 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {561: 5}


class TestSpecialCharacters:
    """Test handling of special characters and delimiters."""

    def test_pipe_delimiter_paired_end(self):
        """Pipe delimiter in paired-end mode."""
        kmer_string = "562:10 |:| 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        # Delimiter should be removed
        assert taxa_dict == {562: 10, 561: 5}
        assert '|:|' not in str(taxa_dict)

    def test_pipe_delimiter_single_end(self):
        """Pipe delimiter when NOT in paired-end mode (shouldn't be removed)."""
        kmer_string = "562:10 |:| 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # When paired_input=False, |:| won't be removed, but should be skipped during parsing
        # because it's malformed (| is not a valid tax_id)
        assert 562 in taxa_dict or 561 in taxa_dict

    def test_reading_frame_delimiter(self):
        """Reading frame delimiters (-:-) should be skipped."""
        kmer_string = "562:10 -:- 561:5 -:- 0:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Reading frame delimiters should be skipped
        assert taxa_dict == {562: 10, 561: 5, 0: 3}
        assert ambiguous_count == 0

    def test_tabs_and_newlines(self):
        """Test that tabs and newlines don't break parsing."""
        kmer_string = "562:10\t561:5\n0:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # split() handles all whitespace
        assert taxa_dict == {562: 10, 561: 5, 0: 3}

    def test_extra_whitespace(self):
        """Test handling of extra whitespace."""
        kmer_string = "562:10    561:5   A:3    "
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10, 561: 5}
        assert ambiguous_count == 3


class TestZeroValues:
    """Test handling of zero and empty values."""

    def test_zero_kmer_count(self):
        """K-mer count of zero."""
        kmer_string = "562:0 561:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 0, 561: 5}
        assert ambiguous_count == 0

    def test_zero_ambiguous_count(self):
        """Ambiguous count of zero."""
        kmer_string = "562:10 A:0"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 10}
        assert ambiguous_count == 0

    def test_all_zero_counts(self):
        """All k-mer counts are zero."""
        kmer_string = "562:0 561:0 A:0"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert sum(taxa_dict.values()) == 0
        assert ambiguous_count == 0

    def test_taxid_zero(self):
        """Tax ID of 0 (unclassified) should be handled."""
        kmer_string = "0:50 A:10"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {0: 50}
        assert ambiguous_count == 10

    def test_empty_string(self):
        """Empty string input."""
        kmer_string = ""
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {}
        assert ambiguous_count == 0

    def test_whitespace_only(self):
        """Only whitespace."""
        kmer_string = "     \t  \n  "
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {}
        assert ambiguous_count == 0


class TestAmbiguousEdgeCases:
    """Edge cases specific to ambiguous k-mer handling (critical for v0.2.0)."""

    def test_ambiguous_only_no_taxa(self):
        """Read with ONLY ambiguous k-mers, no taxa hits."""
        kmer_string = "A:100"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert len(taxa_dict) == 0
        assert ambiguous_count == 100

    def test_ambiguous_with_unclassified_only(self):
        """Ambiguous + unclassified, no real taxa hits."""
        kmer_string = "0:50 A:50"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {0: 50}
        assert ambiguous_count == 50

    def test_many_ambiguous_segments(self):
        """Many separate ambiguous segments."""
        kmer_string = "A:10 562:5 A:5 561:3 A:7 0:2 A:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        assert taxa_dict == {562: 5, 561: 3, 0: 2}
        assert ambiguous_count == 25  # 10+5+7+3

    def test_ambiguous_lowercase_a(self):
        """Test if lowercase 'a' is handled (shouldn't be, but test it)."""
        kmer_string = "562:10 a:5"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=False)

        # Lowercase 'a' should NOT be treated as ambiguous
        # It should fail int('a') and be skipped
        assert taxa_dict == {562: 10}
        assert ambiguous_count == 0


class TestReadClassificationDataclass:
    """Test the ReadClassification dataclass edge cases."""

    def test_read_classification_default_values(self):
        """Test default values of ReadClassification."""
        read = ReadClassification()

        assert read.current_node is None
        assert read.original_conf is None
        assert read.classified is False
        assert read.ambiguous_kmers == 0  # NEW field added in v0.2.0

    def test_read_classification_with_ambiguous_field(self):
        """Test that ambiguous_kmers field exists and works."""
        read = ReadClassification(
            id="test",
            ambiguous_kmers=10
        )

        assert read.ambiguous_kmers == 10

    def test_read_classification_ambiguous_default_zero(self):
        """Test that ambiguous_kmers defaults to 0."""
        read = ReadClassification(id="test")

        assert read.ambiguous_kmers == 0  # Must default to 0, not None


class TestPairedEndEdgeCases:
    """Edge cases specific to paired-end reads."""

    def test_paired_end_empty_read2(self):
        """Paired-end where read2 is empty."""
        kmer_string = "562:50 |:|"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        assert taxa_dict == {562: 50}
        assert ambiguous_count == 0

    def test_paired_end_empty_read1(self):
        """Paired-end where read1 is empty."""
        kmer_string = "|:| 562:50"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        assert taxa_dict == {562: 50}
        assert ambiguous_count == 0

    def test_paired_end_both_empty(self):
        """Paired-end where both reads are empty."""
        kmer_string = "|:|"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        assert taxa_dict == {}
        assert ambiguous_count == 0

    def test_paired_end_multiple_delimiters(self):
        """Multiple |:| delimiters (malformed, but test robustness)."""
        kmer_string = "562:10 |:| 561:5 |:| 0:3"
        taxa_dict, ambiguous_count = process_kmer_string(kmer_string, paired_input=True)

        # remove() only removes first occurrence
        # Second |:| will cause ValueError and be caught
        # But data should still be processed
        assert 562 in taxa_dict or 561 in taxa_dict or 0 in taxa_dict


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
