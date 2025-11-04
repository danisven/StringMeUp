#!/usr/bin/env python3
"""
Pytest fixtures for StringMeUp tests.

Provides reusable test fixtures including mock taxonomy trees,
sample reads, and helper functions.
"""

import pytest
import sys
import os
from dataclasses import dataclass
from typing import Dict, List

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from stringmeup.stringmeup import ReadClassification
from stringmeup.taxonomy import TaxonomyTree, Node


@pytest.fixture
def sample_read_no_ambiguous():
    """
    Sample read without ambiguous k-mers.
    100 k-mers: 60 hits to species 562, 40 unclassified.
    """
    read = ReadClassification(
        id="read1",
        original_taxid=562,
        length="100",
        kmer_string="562:60 0:40",
        classified=False
    )
    return read


@pytest.fixture
def sample_read_with_ambiguous():
    """
    CRITICAL: Sample read WITH ambiguous k-mers.
    100 k-mers: 60 hits to 562, 30 unclassified, 10 ambiguous.

    This is the key test case for the v0.2.0 bug fix.
    """
    read = ReadClassification(
        id="read2",
        original_taxid=562,
        length="100",
        kmer_string="562:60 0:30 A:10",
        classified=False
    )
    return read


@pytest.fixture
def sample_read_high_ambiguous():
    """
    Read with high ambiguous content (50%).
    Tests extreme case where bug had major impact.
    """
    read = ReadClassification(
        id="read3",
        original_taxid=562,
        length="100",
        kmer_string="562:40 0:10 A:50",
        classified=False
    )
    return read


@pytest.fixture
def sample_paired_read_with_ambiguous():
    """
    Paired-end read with ambiguous k-mers in both reads.
    """
    read = ReadClassification(
        id="read4",
        original_taxid=562,
        length="100|95",
        kmer_string="562:50 A:5 |:| 562:30 A:15",
        classified=False
    )
    return read


@pytest.fixture
def sample_read_all_ambiguous():
    """
    Edge case: Read with only ambiguous k-mers.
    Should be marked as unclassified.
    """
    read = ReadClassification(
        id="read5",
        original_taxid=0,
        length="100",
        kmer_string="A:100",
        classified=False
    )
    return read


# Test data for various confidence thresholds
@pytest.fixture(params=[0.0, 0.1, 0.5, 0.6, 0.65, 0.8, 1.0])
def confidence_threshold(request):
    """Parameterized fixture for different confidence thresholds."""
    return request.param


# Test data for ambiguous percentages
@pytest.fixture(params=[0, 5, 10, 20, 30, 50])
def ambiguous_percentage(request):
    """
    Parameterized fixture for different ambiguous k-mer percentages.
    Returns percentage (0-100).
    """
    return request.param


def create_test_read_with_ambiguous_percentage(percent_ambiguous, total_kmers=100, species_taxid=562):
    """
    Helper function to create a test read with specific ambiguous percentage.

    Args:
        percent_ambiguous: Percentage of k-mers that are ambiguous (0-100)
        total_kmers: Total number of k-mers in read
        species_taxid: Tax ID for species hits

    Returns:
        ReadClassification with appropriate k-mer string
    """
    ambiguous_count = int(total_kmers * percent_ambiguous / 100)
    remaining = total_kmers - ambiguous_count

    # Split remaining between species hits (60%) and unclassified (40%)
    species_hits = int(remaining * 0.6)
    unclassified = remaining - species_hits

    kmer_parts = []
    if species_hits > 0:
        kmer_parts.append(f"{species_taxid}:{species_hits}")
    if unclassified > 0:
        kmer_parts.append(f"0:{unclassified}")
    if ambiguous_count > 0:
        kmer_parts.append(f"A:{ambiguous_count}")

    kmer_string = " ".join(kmer_parts)

    return ReadClassification(
        id=f"test_read_{percent_ambiguous}pct_amb",
        original_taxid=species_taxid,
        length=str(total_kmers),
        kmer_string=kmer_string,
        classified=False
    )


# Make helper available to tests
@pytest.fixture
def make_read_with_ambiguous():
    """Fixture that provides the helper function."""
    return create_test_read_with_ambiguous_percentage


# Expected confidence values for regression testing
@pytest.fixture
def expected_confidence_old_vs_new():
    """
    Dictionary mapping ambiguous percentage to (old_confidence, new_confidence).

    These are the expected values for a read with:
    - Total: 100 k-mers
    - Species hits: 60% of non-ambiguous k-mers
    - Unclassified: 40% of non-ambiguous k-mers
    - Ambiguous: varies

    Used for regression testing the v0.2.0 bug fix.
    """
    return {
        0: (0.60, 0.60),    # No change when no ambiguous
        5: (0.632, 0.60),   # Old was inflated by ~3%
        10: (0.667, 0.60),  # Old was inflated by ~7%
        20: (0.75, 0.60),   # Old was inflated by 15%
        30: (0.857, 0.60),  # Old was inflated by 26%
        50: (1.20, 0.60),   # Old was absurdly inflated (>100%!)
    }


@pytest.fixture
def kraken2_example_outputs():
    """
    Real examples of Kraken2 output strings for testing.
    These are actual formats that Kraken2 produces.
    """
    return {
        "single_end_simple": "C\tread1\t562\t100\t562:60 0:40",
        "single_end_with_ambiguous": "C\tread2\t562\t100\t562:60 A:10 0:30",
        "paired_end_simple": "C\tread3\t562\t100|95\t562:50 |:| 562:45",
        "paired_end_with_ambiguous": "C\tread4\t562\t100|95\t562:50 A:5 |:| 562:30 A:10",
        "unclassified": "U\tread5\t0\t100\t0:80 A:20",
        "with_reading_frames": "C\tread6\t562\t100\t562:30 -:- 561:20 -:- 0:50",
    }
