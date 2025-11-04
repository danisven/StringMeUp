#!/usr/bin/env python3
"""
Integration tests for read reclassification using real taxonomy trees.

These tests verify the CORE functionality of StringMeUp:
- Walking up the taxonomy tree
- Calculating confidence at each level
- Stopping when threshold is met
- Handling ambiguous k-mers correctly (v0.2.0 bug fix)

CRITICAL: These tests use a real TaxonomyTree (toy taxonomy) to verify
the complete reclassification workflow. This is what was missing from
the unit tests.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pathlib import Path
from stringmeup.stringmeup import reclassify_read, ReadClassification
from tests.fixtures.toy_taxonomy import create_toy_taxonomy, get_expected_lineage


@pytest.fixture
def toy_taxonomy(tmp_path):
    """
    Create toy taxonomy for testing.

    Taxonomy structure:
    1 (root) → 2 (Bacteria) → 1224 (Proteobacteria) → 561 (Escherichia) → 562 (E. coli)
    """
    return create_toy_taxonomy(tmp_path)


class TestReclassificationBasicBehavior:
    """Test basic reclassification behavior with real taxonomy."""

    def test_low_confidence_moves_up_tree(self, toy_taxonomy):
        """
        Read at E. coli (562) with hits at multiple levels.
        Should move up tree until confidence meets threshold.

        At E. coli: 40/100 = 0.40 < 0.5
        At Escherichia: (40+15)/100 = 0.55 >= 0.5 ✓
        """
        # Create read: 40 hits to E. coli, 15 to Escherichia genus, 45 unclassified
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:40 561:15 0:45"
        )

        # Reclassify with threshold 0.5
        reclassified, lineages = reclassify_read(
            read=read,
            confidence_threshold=0.5,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Should move UP from species (562) to genus (561) or higher
        assert reclassified.reclassified_taxid != 562, \
            "Should move up from E. coli due to low confidence at species level"

        # Should be classified (not unclassified)
        assert reclassified.classified == True
        assert reclassified.reclassified_taxid != 0

        # Confidence at new node should be >= threshold
        assert reclassified.recalculated_conf >= 0.5, \
            f"Final confidence {reclassified.recalculated_conf} should meet threshold 0.5"

        # Should be somewhere in the lineage (likely genus 561)
        lineage = get_expected_lineage(562)
        assert reclassified.reclassified_taxid in lineage, \
            f"Should reclassify to node in lineage: {lineage}"

    def test_high_confidence_stays_at_node(self, toy_taxonomy):
        """
        Read at E. coli with 80% confidence should stay at species level.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:80 0:20"
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

        # Should STAY at E. coli (562)
        assert reclassified.reclassified_taxid == 562, \
            "Should stay at E. coli with 80% confidence"
        assert reclassified.classified == True
        assert reclassified.recalculated_conf >= 0.8

    def test_exact_threshold_passes(self, toy_taxonomy):
        """
        Confidence exactly at threshold should pass (>= check).
        """
        # 50 hits, 50 unclassified = exactly 0.5 confidence
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:50 0:50"
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.5,  # Exactly at threshold
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Should stay at 562 (exactly meets threshold)
        assert reclassified.reclassified_taxid == 562
        assert reclassified.classified == True
        assert reclassified.recalculated_conf == 0.5

    def test_zero_threshold_accepts_anything(self, toy_taxonomy):
        """
        Threshold of 0.0 should accept any classification.
        """
        # Very low confidence
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:10 0:90"
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.0,  # Accept anything
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Should stay at original classification
        assert reclassified.reclassified_taxid == 562
        assert reclassified.classified == True


class TestReclassificationWithAmbiguous:
    """
    CRITICAL TESTS: Verify ambiguous k-mers are included in confidence calculation.
    This is the v0.2.0 bug fix.
    """

    def test_ambiguous_included_in_total(self, toy_taxonomy):
        """
        CRITICAL: Ambiguous k-mers must be included in total for confidence calculation.

        K-mer string: 562:60 0:30 A:10
        Total k-mers: 100 (including 10 ambiguous)
        Confidence at 562: 60/100 = 0.60 (NOT 60/90 = 0.667)
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:60 0:30 A:10"
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.0,  # Accept to check confidence calculation
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Verify ambiguous k-mers were counted
        assert reclassified.ambiguous_kmers == 10, \
            "Must count ambiguous k-mers"

        # Verify total includes ambiguous (60 hits + 30 unclass + 10 ambig = 100)
        # Can't directly check total_kmer_hits, but we can verify confidence
        # If ambiguous were excluded: 60/90 = 0.667
        # If ambiguous were included: 60/100 = 0.60
        assert reclassified.original_conf == 0.60, \
            f"Confidence must be 0.60 (60/100), not {reclassified.original_conf}"

        # Verify it's NOT the buggy value
        buggy_confidence = 60 / 90  # What it would be if ambiguous excluded
        assert reclassified.original_conf < buggy_confidence, \
            "New calculation must be more strict (lower) than buggy version"

    def test_threshold_crossing_due_to_ambiguous(self, toy_taxonomy):
        """
        CRITICAL: Ambiguous k-mers can cause threshold crossing.

        Setup: 60 hits, 30 unclassified, 10 ambiguous, threshold 0.65

        Old buggy: 60/90 = 0.667 → Would PASS threshold → Stay at species
        New correct: 60/100 = 0.60 → FAILS threshold → Move up tree

        This demonstrates the real-world impact of the bug fix.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:60 0:30 A:10"
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.65,  # Above 0.60, below 0.667
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # New correct behavior: should FAIL threshold and move up
        assert reclassified.original_conf == 0.60
        assert reclassified.original_conf < 0.65, "Should be below threshold"

        # Should move up from E. coli
        assert reclassified.reclassified_taxid != 562, \
            "Should move up due to threshold crossing"

        # Old buggy behavior would have stayed at 562
        # We're documenting that this is the CORRECT behavior now

    def test_high_ambiguous_content(self, toy_taxonomy):
        """
        Read with 50% ambiguous k-mers (extreme case).

        Old: 40/50 = 0.80
        New: 40/100 = 0.40

        Massive difference - 40 percentage points!
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:40 0:10 A:50"
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

        assert reclassified.ambiguous_kmers == 50
        assert reclassified.original_conf == 0.40, \
            "With 50% ambiguous, confidence should be 40/100 = 0.40"

        # Old buggy calculation would have been 0.80
        assert reclassified.original_conf < 0.80


class TestReclassificationTreeTraversal:
    """Test that tree traversal works correctly."""

    def test_confidence_increases_up_tree(self, toy_taxonomy):
        """
        As you move up the tree, confidence should increase
        (more k-mers counted in larger clades).

        Read with hits at multiple levels:
        - At E. coli (562): 30 hits
        - At Escherichia (561): +20 hits = 50 total
        - At Proteobacteria (1224): +10 hits = 60 total
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,  # Start at E. coli
            kmer_string="562:30 561:20 1224:10 0:40"
        )

        # At E. coli: 30/100 = 0.30
        # At Escherichia: (30+20)/100 = 0.50
        # At Proteobacteria: (30+20+10)/100 = 0.60

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.55,  # Between 0.50 and 0.60
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Should reclassify to Proteobacteria (1224) or higher
        assert reclassified.reclassified_taxid in [1224, 2, 1], \
            "Should move up to Proteobacteria or higher"
        assert reclassified.recalculated_conf >= 0.55

    def test_stops_at_root(self, toy_taxonomy):
        """
        If confidence never meets threshold, should stop at root (1).
        """
        # Only 20% hits - will never meet 0.5 threshold
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:20 0:80"
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

        # Should reach root and be unclassified
        assert reclassified.classified == False
        assert reclassified.reclassified_taxid == 0


class TestReclassificationMaxConfidence:
    """Test max_confidence calculation and early termination."""

    def test_max_confidence_calculated(self, toy_taxonomy):
        """
        Max confidence should be total hits / total kmers.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:60 0:40"  # 60 hits, 40 unclassified
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

        # Max confidence should be 60/100 = 0.60
        assert reclassified.max_confidence == 0.60

    def test_doomed_to_fail_early_exit(self, toy_taxonomy):
        """
        If max confidence < threshold, should exit early without tree walking.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:30 0:70"  # Max 30% confidence
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.5,  # Impossible to meet
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=False
        )

        # Should be unclassified (doomed to fail)
        assert reclassified.classified == False
        assert reclassified.reclassified_taxid == 0
        assert reclassified.max_confidence == 0.30

        # Should have set max_confidence
        assert reclassified.recalculated_conf == 0.30


class TestReclassificationPairedEnd:
    """Test reclassification with paired-end reads."""

    def test_paired_end_with_ambiguous_both_reads(self, toy_taxonomy):
        """
        Paired-end: ambiguous k-mers from both reads must be aggregated.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:40 A:5 |:| 562:30 A:5"  # 70 hits, 10 ambiguous
        )

        reclassified, _ = reclassify_read(
            read=read,
            confidence_threshold=0.0,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages={},
            paired_input=True  # PAIRED
        )

        # Should aggregate both reads
        assert reclassified.ambiguous_kmers == 10, "Ambiguous from both reads: 5+5"
        # Total: 70 hits + 10 ambiguous = 80
        # Confidence: 70/80 = 0.875
        assert reclassified.original_conf == 0.875


class TestReclassificationLineageCaching:
    """Test that lineage caching works correctly."""

    def test_lineage_cache_populated(self, toy_taxonomy):
        """
        Lineage cache should be populated during reclassification.
        """
        read = ReadClassification(
            id="test_read",
            original_taxid=562,
            kmer_string="562:60 0:40"
        )

        lineage_cache = {}
        reclassified, new_cache = reclassify_read(
            read=read,
            confidence_threshold=0.0,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages=lineage_cache,
            paired_input=False
        )

        # Cache should have been populated
        assert len(new_cache) > 0
        assert 562 in new_cache  # E. coli lineage should be cached

    def test_lineage_cache_reused(self, toy_taxonomy):
        """
        Lineage cache should be reused on subsequent calls.
        """
        # First call - populate cache
        read1 = ReadClassification(id="read1", original_taxid=562, kmer_string="562:60 0:40")
        lineage_cache = {}

        _, cache1 = reclassify_read(
            read=read1,
            confidence_threshold=0.0,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages=lineage_cache,
            paired_input=False
        )

        # Second call - reuse cache
        read2 = ReadClassification(id="read2", original_taxid=561, kmer_string="561:50 0:50")

        _, cache2 = reclassify_read(
            read=read2,
            confidence_threshold=0.0,
            taxonomy_tree=toy_taxonomy,
            verbose_input=False,
            minimum_hit_groups=1,
            taxa_lineages=cache1,  # Reuse
            paired_input=False
        )

        # Cache should have entries from both calls
        assert 562 in cache2
        assert 561 in cache2


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
