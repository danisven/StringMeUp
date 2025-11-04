#!/usr/bin/env python3
"""
Tests for dataclasses and data structures.

Tests the dataclass definitions and their behavior.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from stringmeup.stringmeup import ReadClassification, ReportNode


class TestReadClassificationDataclass:
    """Comprehensive tests for ReadClassification dataclass."""

    def test_default_initialization(self):
        """Test ReadClassification with default values."""
        read = ReadClassification()

        # Check all default values
        assert read.id is None
        assert read.classified == False
        assert read.original_taxid is None
        assert read.reclassified_taxid is None
        assert read.current_node is None
        assert read.kmer_string is None
        assert read.original_conf is None
        assert read.recalculated_conf is None
        assert read.max_confidence is None
        assert read.minimizer_hit_groups is None
        assert read.ambiguous_kmers == 0  # CRITICAL: v0.2.0 field

    def test_partial_initialization(self):
        """Test ReadClassification with some fields set."""
        read = ReadClassification(
            id="read1",
            classified=True,
            original_taxid=562,
            kmer_string="562:100"
        )

        assert read.id == "read1"
        assert read.classified == True
        assert read.original_taxid == 562
        assert read.kmer_string == "562:100"
        # Others should still be defaults
        assert read.reclassified_taxid is None
        assert read.ambiguous_kmers == 0

    def test_full_initialization(self):
        """Test ReadClassification with all major fields set."""
        read = ReadClassification(
            id="read2",
            classified=True,
            original_taxid=562,
            reclassified_taxid=2,
            current_node=2,
            kmer_string="562:60 0:30 A:10",
            original_conf=0.60,
            recalculated_conf=0.80,
            ambiguous_kmers=10
        )

        # Set non-typed fields after initialization
        read.max_confidence = 0.60
        read.minimizer_hit_groups = 5

        assert read.id == "read2"
        assert read.classified == True
        assert read.original_taxid == 562
        assert read.reclassified_taxid == 2
        assert read.current_node == 2
        assert read.kmer_string == "562:60 0:30 A:10"
        assert read.original_conf == 0.60
        assert read.recalculated_conf == 0.80
        assert read.max_confidence == 0.60
        assert read.minimizer_hit_groups == 5
        assert read.ambiguous_kmers == 10

    def test_ambiguous_kmers_field_exists(self):
        """
        CRITICAL v0.2.0 TEST: Verify ambiguous_kmers field exists.
        This field was added in v0.2.0 to fix the confidence bug.
        """
        read = ReadClassification()

        # Field should exist and default to 0
        assert hasattr(read, 'ambiguous_kmers')
        assert read.ambiguous_kmers == 0

    def test_field_mutation(self):
        """Test that fields can be mutated after initialization."""
        read = ReadClassification(id="read3")

        # Mutate fields
        read.classified = True
        read.original_taxid = 562
        read.ambiguous_kmers = 15

        assert read.classified == True
        assert read.original_taxid == 562
        assert read.ambiguous_kmers == 15

    def test_classification_workflow(self):
        """
        Test a typical classification workflow using the dataclass.
        Simulates the workflow from unclassified to classified.
        """
        # Start with unclassified read
        read = ReadClassification(
            id="read4",
            kmer_string="562:60 0:30 A:10",
            original_taxid=562,
            classified=False
        )

        # Simulate processing
        read.ambiguous_kmers = 10
        read.original_conf = 0.60
        read.max_confidence = 0.60

        # Simulate reclassification to parent
        read.current_node = 2  # Moved to parent
        read.recalculated_conf = 0.80
        read.classified = True
        read.reclassified_taxid = 2

        # Verify final state
        assert read.classified == True
        assert read.reclassified_taxid == 2
        assert read.recalculated_conf > read.original_conf
        assert read.ambiguous_kmers == 10

    def test_confidence_values_none_vs_zero(self):
        """Test that confidence values can be None (unset) vs 0.0 (set to zero)."""
        read1 = ReadClassification()
        assert read1.original_conf is None  # Not yet set

        read2 = ReadClassification(original_conf=0.0)
        assert read2.original_conf == 0.0  # Explicitly set to zero

        # These should be different
        assert read1.original_conf is not read2.original_conf


class TestReportNodeDataclass:
    """Tests for ReportNode dataclass."""

    def test_minimal_initialization(self):
        """Test ReportNode with minimal required fields."""
        # ReportNode requires all fields (no defaults)
        node = ReportNode(
            ratio="0.00",
            hits_at_clade=0,
            hits_at_node=0,
            rank_code="U",
            rank_depth=0,
            node_taxid=0,
            name="test",
            offset=0
        )

        assert node.ratio == "0.00"
        assert node.hits_at_clade == 0
        assert node.hits_at_node == 0
        assert node.rank_code == "U"
        assert node.rank_depth == 0
        assert node.node_taxid == 0
        assert node.name == "test"
        assert node.offset == 0

    def test_full_initialization(self):
        """Test ReportNode with all fields set."""
        node = ReportNode(
            ratio="25.50",
            hits_at_clade=255,
            hits_at_node=100,
            rank_code="S",
            rank_depth=7,
            node_taxid=562,
            offset=2,
            name="Escherichia coli"
        )

        assert node.ratio == "25.50"
        assert node.hits_at_clade == 255
        assert node.hits_at_node == 100
        assert node.rank_code == "S"
        assert node.rank_depth == 7
        assert node.node_taxid == 562
        assert node.offset == 2
        assert node.name == "Escherichia coli"

    def test_unclassified_report_node(self):
        """Test ReportNode for unclassified reads (tax_id 0)."""
        node = ReportNode(
            ratio="10.00",
            hits_at_clade=100,
            hits_at_node=100,
            rank_code="U",
            rank_depth=0,
            node_taxid=0,
            offset=0,
            name="unclassified"
        )

        assert node.rank_code == "U"
        assert node.node_taxid == 0
        assert node.name == "unclassified"

    def test_root_report_node(self):
        """Test ReportNode for root (tax_id 1)."""
        node = ReportNode(
            ratio="90.00",
            hits_at_clade=900,
            hits_at_node=0,
            rank_code="R",
            rank_depth=0,
            node_taxid=1,
            offset=0,
            name="root"
        )

        assert node.node_taxid == 1
        assert node.name == "root"
        assert node.hits_at_node == 0  # Root usually has 0 direct hits

    def test_species_level_report_node(self):
        """Test ReportNode for species-level classification."""
        node = ReportNode(
            ratio="5.25",
            hits_at_clade=52,
            hits_at_node=52,
            rank_code="S",
            rank_depth=7,
            node_taxid=562,
            offset=6,
            name="Escherichia coli"
        )

        assert node.rank_code == "S"
        assert node.hits_at_clade == node.hits_at_node  # Leaf node


class TestDataclassFieldTypes:
    """Test that dataclass fields accept correct types."""

    def test_read_classification_string_fields(self):
        """Test string field types."""
        read = ReadClassification(
            id="test_id",
            kmer_string="562:10 A:5"
        )

        assert isinstance(read.id, str)
        assert isinstance(read.kmer_string, str)

    def test_read_classification_int_fields(self):
        """Test integer field types."""
        read = ReadClassification(
            original_taxid=562,
            reclassified_taxid=2,
            current_node=2,
            ambiguous_kmers=10
        )

        # Set non-typed field
        read.minimizer_hit_groups = 5

        assert isinstance(read.original_taxid, int)
        assert isinstance(read.reclassified_taxid, int)
        assert isinstance(read.current_node, int)
        assert isinstance(read.minimizer_hit_groups, int)
        assert isinstance(read.ambiguous_kmers, int)

    def test_read_classification_float_fields(self):
        """Test float field types."""
        read = ReadClassification(
            original_conf=0.60,
            recalculated_conf=0.80
        )

        # Set non-typed field
        read.max_confidence = 0.60

        assert isinstance(read.original_conf, float)
        assert isinstance(read.recalculated_conf, float)
        assert isinstance(read.max_confidence, float)

    def test_read_classification_bool_field(self):
        """Test boolean field type."""
        read = ReadClassification(classified=True)

        assert isinstance(read.classified, bool)
        assert read.classified == True


class TestDataclassRegressionTests:
    """Regression tests for dataclass behavior related to v0.2.0 bug fix."""

    def test_ambiguous_kmers_zero_not_none(self):
        """
        CRITICAL: ambiguous_kmers must default to 0, not None.
        This is important for confidence calculations.
        """
        read = ReadClassification()

        # Must be 0, not None
        assert read.ambiguous_kmers == 0
        assert read.ambiguous_kmers is not None

        # Should work in arithmetic without None checks
        total = 100
        confidence = 60 / (total + read.ambiguous_kmers)
        assert confidence == 0.60  # Should not raise TypeError

    def test_confidence_calculation_with_ambiguous(self):
        """
        Test that confidence can be calculated correctly using ambiguous_kmers.
        This is the core of the v0.2.0 bug fix.
        """
        read = ReadClassification(
            kmer_string="562:60 0:30 A:10",
            ambiguous_kmers=10
        )

        # Simulate calculating total from kmer string
        taxa_kmers = 60 + 30  # From "562:60 0:30"
        total_kmers = taxa_kmers + read.ambiguous_kmers  # +10 from A:10

        # Calculate confidence (hits / total including ambiguous)
        hits = 60
        confidence = hits / total_kmers

        assert total_kmers == 100
        assert confidence == 0.60
        assert read.ambiguous_kmers == 10


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
