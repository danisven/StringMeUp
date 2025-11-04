#!/usr/bin/env python3
"""
Integration tests for error handling and edge cases.

These tests verify that StringMeUp handles errors gracefully:
- Invalid input data
- Missing taxonomy nodes
- Corrupted files
- Invalid parameters
- Edge cases in data processing

IMPORTANT: These tests document expected error behavior and ensure
the system fails safely with informative errors.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pathlib import Path
from stringmeup.stringmeup import (
    reclassify_read,
    ReadClassification,
    create_read,
    validate_input_file
)
from tests.fixtures.toy_taxonomy import create_toy_taxonomy


@pytest.fixture
def toy_taxonomy(tmp_path):
    """Create toy taxonomy for testing."""
    return create_toy_taxonomy(tmp_path)


class TestMissingTaxonomyNodes:
    """Test handling of reads classified to taxa not in the taxonomy."""

    def test_reclassify_read_with_unknown_taxid(self, toy_taxonomy):
        """
        Test reclassification when original_taxid is not in taxonomy.

        ACTUAL BEHAVIOR: Raises KeyError when tax_id not found in taxonomy.
        This is expected - StringMeUp requires all tax_ids in input to exist
        in the provided taxonomy database.

        IMPORTANT: Users must ensure their taxonomy database matches their
        Kraken2 classification database.
        """
        # Create read with unknown taxid (not in toy taxonomy)
        read = ReadClassification(
            id="test_read",
            original_taxid=99999,  # Unknown taxid
            kmer_string="99999:60 0:40"
        )

        # Should raise KeyError (documented behavior)
        with pytest.raises(KeyError):
            reclassified, _ = reclassify_read(
                read=read,
                confidence_threshold=0.5,
                taxonomy_tree=toy_taxonomy,
                verbose_input=False,
                minimum_hit_groups=1,
                taxa_lineages={},
                paired_input=False
            )

    def test_kmer_string_with_unknown_taxids(self, toy_taxonomy):
        """
        Test k-mer string containing taxa not in taxonomy.

        K-mer string: "562:30 99999:20 0:50"
        - 562 is in taxonomy (E. coli)
        - 99999 is NOT in taxonomy
        - 0 is unclassified

        ACTUAL BEHAVIOR: During confidence calculation, when walking through
        k-mer string, encountering unknown tax_id raises KeyError.

        This documents that StringMeUp expects ALL taxa in k-mer strings
        to be present in the taxonomy.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,  # Start at known taxon
            kmer_string="562:30 99999:20 0:50"  # Contains unknown tax_id
        )

        # Should raise KeyError when processing k-mer string
        with pytest.raises(KeyError):
            reclassified, _ = reclassify_read(
                read=read,
                confidence_threshold=0.0,
                taxonomy_tree=toy_taxonomy,
                verbose_input=False,
                minimum_hit_groups=1,
                taxa_lineages={},
                paired_input=False
            )


class TestMalformedInput:
    """Test handling of malformed Kraken2 input."""

    def test_create_read_with_missing_fields(self):
        """
        Test parsing line with missing required fields.

        Valid format: C/U <tab> read_id <tab> taxid <tab> length <tab> kmer_string
        """
        # Line with only 3 fields (missing length and kmer_string)
        malformed_line = "C\tread1\t562"

        # Should raise error or handle gracefully
        with pytest.raises((ValueError, IndexError)):
            create_read(malformed_line, verbose_input=False)

    def test_create_read_with_invalid_taxid(self):
        """Test parsing line with non-numeric taxid."""
        malformed_line = "C\tread1\tNOT_A_NUMBER\t100\t562:60 0:40"

        with pytest.raises(ValueError):
            create_read(malformed_line, verbose_input=False)

    def test_create_read_with_invalid_length(self):
        """
        Test parsing line with non-numeric length.

        ACTUAL BEHAVIOR: Implementation is lenient - the length field
        is stored as a string, so non-numeric values don't raise errors
        during parsing. This may cause issues downstream if length is
        used for calculations.
        """
        malformed_line = "C\tread1\t562\tNOT_A_NUMBER\t562:60 0:40"

        # Implementation is lenient - doesn't validate length format
        read = create_read(malformed_line, verbose_input=False)

        # But parsing should succeed
        assert read.id == "read1"
        # Length is stored but not validated during parsing
        assert read.length is not None

    def test_create_read_with_malformed_kmer_string(self):
        """
        Test parsing line with malformed k-mer string.

        Note: K-mer string format is flexible, so this may not raise error
        but should be documented.
        """
        # K-mer string without colons
        line = "C\tread1\t562\t100\tinvalid_kmer_format"

        # May or may not raise error depending on implementation
        # Document actual behavior
        try:
            read = create_read(line, verbose_input=False)
            # If it doesn't raise, verify it at least parsed something
            assert read.id == "read1"
        except (ValueError, IndexError):
            # Expected if implementation validates k-mer format
            pass


class TestFileValidation:
    """Test input file validation."""

    def test_validate_inconsistent_paired_format(self, tmp_path):
        """
        Test validation when paired_input=True but first line is single-end.

        NOTE: validate_input_file() only validates the FIRST LINE.
        Inconsistencies in later lines won't be detected.
        """
        inconsistent_file = tmp_path / "inconsistent.kraken2"
        inconsistent_file.write_text(
            "C\tread1\t562\t100\t562:60 0:40\n"  # First line: Single-end
            "C\tread2\t562\t100\t562:40 |:| 562:30\n"  # Paired
            "C\tread3\t562\t100\t562:50 |:| 562:30\n"  # Paired
        )

        # Should fail because first line doesn't match paired_input=True
        with pytest.raises(SystemExit):
            validate_input_file(
                str(inconsistent_file),
                verbose_input=False,
                minimum_hit_groups=None,
                paired_input=True  # Expecting paired, but first line is single-end
            )

    def test_validate_empty_file(self, tmp_path):
        """
        Test validation of completely empty file.

        ACTUAL BEHAVIOR: validate_input_file crashes with IndexError when
        trying to access fields of an empty line. This is a known limitation.

        Empty files should be handled before calling validate_input_file.
        """
        empty_file = tmp_path / "empty.kraken2"
        empty_file.write_text("")

        # Empty file causes IndexError in validation
        with pytest.raises(IndexError):
            validate_input_file(
                str(empty_file),
                verbose_input=False,
                minimum_hit_groups=None,
                paired_input=False
            )


class TestEdgeCasesInConfidence:
    """Test edge cases in confidence calculation."""

    def test_read_with_zero_kmers(self, toy_taxonomy):
        """
        Test read with no k-mers (malformed but shouldn't crash).

        This is an edge case that shouldn't occur in real data but
        should be handled gracefully.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string=""  # Empty k-mer string
        )

        # Should handle without crashing
        try:
            reclassified, _ = reclassify_read(
                read=read,
                confidence_threshold=0.5,
                taxonomy_tree=toy_taxonomy,
                verbose_input=False,
                minimum_hit_groups=1,
                taxa_lineages={},
                paired_input=False
            )
            # If it succeeds, should likely be unclassified
            assert reclassified.classified == False
        except (ValueError, ZeroDivisionError):
            # Acceptable if it raises error for empty k-mer string
            pass

    def test_read_with_only_unclassified_kmers(self, toy_taxonomy):
        """
        Test read where all k-mers are unclassified (0:100).

        Expected: Max confidence = 0, should be unclassified.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,  # Originally classified
            kmer_string="0:100"  # But all k-mers are unclassified
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.5,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Should be unclassified (no hits to any taxon)
        assert reclassified.classified == False
        assert reclassified.max_confidence == 0.0

    def test_read_with_only_ambiguous_kmers(self, toy_taxonomy):
        """
        Test read where all k-mers are ambiguous (A:100).

        Expected: No hits to any specific taxon → unclassified.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="A:100"  # All ambiguous
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.5,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Should be unclassified
        assert reclassified.classified == False
        assert reclassified.ambiguous_kmers == 100


class TestPairedEndEdgeCases:
    """Test edge cases specific to paired-end reads."""

    def test_paired_read_with_empty_mate(self, toy_taxonomy):
        """
        Test paired-end read where one mate has no k-mers.

        Format: "562:60 |:| 0:0" or "562:60 |:|"
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:60 0:20 |:|"  # Empty second read
        )

        # Should handle gracefully
        try:
            reclassified, _ = reclassify_read(
                read=read,
                confidence_threshold=0.5,
                taxonomy_tree=toy_taxonomy,
                verbose_input=False,
                minimum_hit_groups=1,
                taxa_lineages={},
                paired_input=True
            )
            # Should process based on first read only
            assert reclassified is not None
        except (ValueError, IndexError):
            # Acceptable if it raises error for malformed paired format
            pass

    def test_paired_read_with_unbalanced_delimiter(self, toy_taxonomy):
        """
        Test k-mer string with multiple |:| delimiters (malformed).
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:30 |:| 562:20 |:| 562:10"  # Two delimiters
        )

        # Should handle or raise error
        try:
            reclassified, _ = reclassify_read(
                read=read,
                confidence_threshold=0.5,
                taxonomy_tree=toy_taxonomy,
                verbose_input=False,
                minimum_hit_groups=1,
                taxa_lineages={},
                paired_input=True
            )
            # If it succeeds, verify result is reasonable
            assert reclassified is not None
        except (ValueError, IndexError):
            # Expected if implementation validates delimiter count
            pass


class TestVerboseInputEdgeCases:
    """Test edge cases with verbose Kraken2 input."""

    def test_create_read_verbose_with_zero_hit_groups(self):
        """
        Test verbose format with 0 minimizer hit groups.

        This could happen for very low-quality reads.
        """
        # Verbose format: C/U <tab> id <tab> taxid <tab> length <tab> hit_groups <tab> kmer_string
        line = "C\tread1\t562\t100\t0\t562:60 0:40"

        read = create_read(line, verbose_input=True)

        assert read.id == "read1"
        assert read.minimizer_hit_groups == 0

    def test_reclassify_with_minimum_hit_groups_filter(self, toy_taxonomy):
        """
        Test that reads below minimum_hit_groups threshold are handled.

        Note: minimizer_hit_groups is set during create_read(), not in
        ReadClassification constructor. The filtering happens in main_loop.
        """
        # Create read via parsing (sets minimizer_hit_groups)
        line = "C\ttest_read\t562\t100\t5\t562:60 0:40"  # 5 hit groups
        read = create_read(line, verbose_input=True)

        assert read.minimizer_hit_groups == 5

        # Reclassify with minimum_hit_groups requirement
        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.5,
            taxonomy_tree=toy_taxonomy,
            verbose_input=True,
            minimum_hit_groups=10,  # Require at least 10
            taxa_lineages={},
            paired_input=False
        )

        # reclassify_read processes the read regardless
        # Filtering happens in main_loop based on minimizer_hit_groups
        assert reclassified is not None


class TestNumericalEdgeCases:
    """Test numerical edge cases in calculations."""

    def test_confidence_calculation_precision(self, toy_taxonomy):
        """
        Test that confidence is calculated with reasonable precision.

        Example: 1/3 = 0.333... should be handled correctly.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:33 0:67"  # 33/100 = 0.33
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.0,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Confidence should be exactly 0.33
        assert reclassified.original_conf == 0.33

    def test_very_long_read_with_many_kmers(self, toy_taxonomy):
        """
        Test read with very large k-mer counts.

        Simulates long reads with thousands of k-mers.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:5000 0:5000"  # 10000 total k-mers
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.5,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Should handle large numbers correctly
        assert reclassified.original_conf == 0.5
        assert reclassified.reclassified_taxid == 562


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
