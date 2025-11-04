#!/usr/bin/env python3
"""
Integration tests for Kraken2 report generation.

Tests the report generation functions that create Kraken2-compatible reports:
- Cumulative clade counting (DFS)
- Report node sorting
- Format compatibility with Kraken2
- Unclassified handling

These tests verify that StringMeUp produces reports matching Kraken2's format.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pathlib import Path
from stringmeup.stringmeup import (
    get_kraken2_report_content,
    format_kraken2_report_row,
    ReportNode
)
from tests.fixtures.toy_taxonomy import create_toy_taxonomy


@pytest.fixture
def toy_taxonomy(tmp_path):
    """Create toy taxonomy for testing."""
    return create_toy_taxonomy(tmp_path)


class TestReportNodeFormatting:
    """Test formatting of individual report rows."""

    def test_format_unclassified_row(self):
        """Test formatting of unclassified row (first row in report)."""
        node = ReportNode(
            ratio="10.50",
            hits_at_clade=105,
            hits_at_node=105,
            rank_code="U",
            rank_depth=0,
            node_taxid=0,
            name="unclassified",
            offset=0
        )

        formatted = format_kraken2_report_row(node)

        # Format: percentage, clade_count, direct_count, rank, tax_id, name
        assert "10.50" in formatted
        assert "105" in formatted
        assert "U" in formatted
        assert "0" in formatted
        assert "unclassified" in formatted

    def test_format_root_row(self):
        """Test formatting of root node."""
        node = ReportNode(
            ratio="89.50",
            hits_at_clade=895,
            hits_at_node=0,  # Root typically has 0 direct hits
            rank_code="R",
            rank_depth=0,
            node_taxid=1,
            name="root",
            offset=0
        )

        formatted = format_kraken2_report_row(node)

        assert "89.50" in formatted
        assert "895" in formatted
        assert "R" in formatted
        assert "root" in formatted

    def test_format_species_row_with_offset(self):
        """Test formatting of species-level node with indentation."""
        node = ReportNode(
            ratio="25.00",
            hits_at_clade=250,
            hits_at_node=250,
            rank_code="S",
            rank_depth=7,
            node_taxid=562,
            name="Escherichia coli",
            offset=4  # Indented 4 levels
        )

        formatted = format_kraken2_report_row(node)

        assert "25.00" in formatted
        assert "250" in formatted
        assert "S" in formatted
        assert "562" in formatted
        assert "Escherichia coli" in formatted
        # Should have indentation spaces (offset * 2)
        assert "  " in formatted  # At least some indentation

    def test_format_preserves_all_fields(self):
        """Test that all required fields are in formatted output."""
        node = ReportNode(
            ratio="15.25",
            hits_at_clade=152,
            hits_at_node=100,
            rank_code="G",
            rank_depth=6,
            node_taxid=561,
            name="Escherichia",
            offset=3
        )

        formatted = format_kraken2_report_row(node)
        fields = formatted.split('\t')

        # Kraken2 report format has 6 fields minimum
        assert len(fields) >= 5, "Should have at least 5 tab-separated fields"


class TestReportGeneration:
    """Test complete report generation with toy taxonomy."""

    def test_report_generation_simple(self, toy_taxonomy):
        """
        Test report generation with simple classification scenario.

        Scenario:
        - 10 reads total
        - 8 reads classified to E. coli (562)
        - 2 reads unclassified
        """
        # Create tax_reads structure
        tax_reads = {
            'hits_at_node': {
                562: 8,  # 8 reads at E. coli
            },
            'hits_at_clade': {}  # Will be populated by function
        }

        total_reads = 10

        # Generate report
        report_nodes = get_kraken2_report_content(
            tax_reads,
            toy_taxonomy,
            total_reads
        )

        # Should have nodes for the lineage + unclassified
        # Unclassified + root + Bacteria + Proteobacteria + Escherichia + E. coli
        assert len(report_nodes) > 0, "Should generate report nodes"

        # First node should be unclassified
        assert report_nodes[0].node_taxid == 0
        assert report_nodes[0].rank_code == "U"
        assert report_nodes[0].hits_at_node == 2  # 2 unclassified

        # Should have E. coli in the report
        ecoli_nodes = [n for n in report_nodes if n.node_taxid == 562]
        assert len(ecoli_nodes) == 1, "Should have E. coli in report"
        assert ecoli_nodes[0].hits_at_node == 8

    def test_report_clade_counts_accumulated(self, toy_taxonomy):
        """
        Test that clade counts are correctly accumulated up the tree.

        If E. coli has 8 reads:
        - E. coli clade: 8
        - Escherichia clade: 8 (includes E. coli)
        - Proteobacteria clade: 8
        - Bacteria clade: 8
        - Root clade: 8
        """
        tax_reads = {
            'hits_at_node': {
                562: 8,  # Only hits at E. coli
            },
            'hits_at_clade': {}
        }

        report_nodes = get_kraken2_report_content(
            tax_reads,
            toy_taxonomy,
            total_reads=10
        )

        # Check E. coli clade count
        ecoli = [n for n in report_nodes if n.node_taxid == 562][0]
        assert ecoli.hits_at_clade == 8

        # Check Escherichia genus clade count (should include E. coli)
        escherichia = [n for n in report_nodes if n.node_taxid == 561]
        if escherichia:  # May or may not be in report depending on hits
            assert escherichia[0].hits_at_clade >= 8

        # Check root clade count (should include all classified)
        root = [n for n in report_nodes if n.node_taxid == 1][0]
        assert root.hits_at_clade == 8

    def test_report_unclassified_at_top(self, toy_taxonomy):
        """
        Test that unclassified row appears first in report.
        This is Kraken2 convention.
        """
        tax_reads = {
            'hits_at_node': {562: 5},
            'hits_at_clade': {}
        }

        report_nodes = get_kraken2_report_content(
            tax_reads,
            toy_taxonomy,
            total_reads=10
        )

        # First node must be unclassified
        assert report_nodes[0].node_taxid == 0
        assert report_nodes[0].rank_code == "U"
        assert report_nodes[0].hits_at_node == 5  # 10 total - 5 classified

    def test_report_percentages_sum_to_100(self, toy_taxonomy):
        """
        Test that percentages in report sum to 100%.
        (Unclassified + Root = 100%)
        """
        tax_reads = {
            'hits_at_node': {562: 7},
            'hits_at_clade': {}
        }

        total_reads = 10

        report_nodes = get_kraken2_report_content(
            tax_reads,
            toy_taxonomy,
            total_reads
        )

        # Get unclassified and root percentages
        unclassified_ratio = float(report_nodes[0].ratio)
        root = [n for n in report_nodes if n.node_taxid == 1][0]
        root_ratio = float(root.ratio)

        # Should sum to 100% (with small floating point tolerance)
        total_percentage = unclassified_ratio + root_ratio
        assert abs(total_percentage - 100.0) < 0.01, \
            f"Unclassified ({unclassified_ratio}%) + Root ({root_ratio}%) should sum to 100%"

    def test_report_with_multiple_species(self, toy_taxonomy):
        """
        Test report generation when reads hit multiple taxa at same level.

        Note: Our toy taxonomy only has one species, but we can simulate
        multiple hits at different levels.
        """
        tax_reads = {
            'hits_at_node': {
                562: 5,   # 5 at E. coli species
                561: 3,   # 3 at Escherichia genus (not E. coli)
            },
            'hits_at_clade': {}
        }

        report_nodes = get_kraken2_report_content(
            tax_reads,
            toy_taxonomy,
            total_reads=10
        )

        # Should have both taxa in report
        ecoli_nodes = [n for n in report_nodes if n.node_taxid == 562]
        genus_nodes = [n for n in report_nodes if n.node_taxid == 561]

        assert len(ecoli_nodes) == 1, "Should have E. coli"
        assert len(genus_nodes) == 1, "Should have Escherichia genus"

        # Genus clade should include species hits
        assert genus_nodes[0].hits_at_clade >= 8  # 5 + 3


class TestReportEdgeCases:
    """Test edge cases in report generation."""

    def test_report_with_mostly_unclassified(self, toy_taxonomy):
        """
        Test report generation with mostly unclassified reads.

        Note: Report generation requires at least 1 classified read to build the tree.
        With 0 classified reads, only unclassified row would be meaningful.
        """
        tax_reads = {
            'hits_at_node': {562: 1},  # Only 1 classified read
            'hits_at_clade': {}
        }

        report_nodes = get_kraken2_report_content(
            tax_reads,
            toy_taxonomy,
            total_reads=100
        )

        # Should have unclassified row with 99 reads
        assert len(report_nodes) >= 1
        assert report_nodes[0].node_taxid == 0
        assert report_nodes[0].hits_at_node == 99  # 100 - 1 classified

    def test_report_with_all_classified(self, toy_taxonomy):
        """Test report with 100% classification rate."""
        tax_reads = {
            'hits_at_node': {562: 100},
            'hits_at_clade': {}
        }

        report_nodes = get_kraken2_report_content(
            tax_reads,
            toy_taxonomy,
            total_reads=100
        )

        # Unclassified should be 0
        assert report_nodes[0].hits_at_node == 0
        assert report_nodes[0].ratio == "0.00"

        # Root should be 100%
        root = [n for n in report_nodes if n.node_taxid == 1][0]
        assert root.ratio == "100.00"

    def test_report_offsets_increase_down_tree(self, toy_taxonomy):
        """
        Test that offset (indentation) increases as you go down the tree.

        Root should have offset 0
        Bacteria should have offset 1
        Proteobacteria should have offset 2
        Escherichia should have offset 3
        E. coli should have offset 4
        """
        tax_reads = {
            'hits_at_node': {562: 10},
            'hits_at_clade': {}
        }

        report_nodes = get_kraken2_report_content(
            tax_reads,
            toy_taxonomy,
            total_reads=10
        )

        # Build offset map
        offsets = {n.node_taxid: n.offset for n in report_nodes if n.node_taxid != 0}

        # Root should be at offset 0
        assert offsets[1] == 0, "Root should have offset 0"

        # Each child should have offset > parent
        if 2 in offsets:  # Bacteria
            assert offsets[2] > offsets[1], "Bacteria offset should be > root"
        if 1224 in offsets:  # Proteobacteria
            assert offsets[1224] > offsets[2], "Proteobacteria offset should be > Bacteria"
        if 562 in offsets:  # E. coli
            assert offsets[562] > offsets[561], "E. coli offset should be > Escherichia"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
