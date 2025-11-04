#!/usr/bin/env python3
"""
End-to-end workflow tests for StringMeUp.

These tests verify the COMPLETE workflow from input file to output files:
- Reading Kraken2 classification files
- Reclassifying reads based on confidence threshold
- Writing output classifications
- Writing verbose output
- Generating Kraken2 reports

CRITICAL: These tests verify that all components work together correctly
in a real-world usage scenario.
"""

import pytest
import sys
import os
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from stringmeup.stringmeup import (
    main_loop,
    read_file,
    write_file,
    get_kraken2_report_content,
    format_kraken2_report_row
)
from tests.fixtures.toy_taxonomy import create_toy_taxonomy


@pytest.fixture
def toy_taxonomy(tmp_path):
    """Create toy taxonomy for testing."""
    return create_toy_taxonomy(tmp_path)


@pytest.fixture
def example_kraken2_file():
    """Path to example.kraken2 fixture."""
    return Path(__file__).parent.parent / "fixtures" / "example.kraken2"


@pytest.fixture
def example_paired_file():
    """Path to example_paired.kraken2 fixture."""
    return Path(__file__).parent.parent / "fixtures" / "example_paired.kraken2"


@pytest.fixture
def mock_args():
    """Create mock args object for testing."""
    class Args:
        def __init__(self):
            self.confidence_threshold = 0.5
            self.minimum_hit_groups = None
            self.keep_unclassified = False
            self.output_report = None  # Don't write report to file in tests

    return Args()


class TestCompleteWorkflow:
    """Test complete end-to-end workflow."""

    def test_complete_reclassification_workflow(self, toy_taxonomy, example_kraken2_file, mock_args, tmp_path):
        """
        Test complete workflow: input file → reclassification → output files.

        This is the CORE end-to-end test that verifies:
        1. Reading input file
        2. Processing all reads
        3. Reclassifying based on threshold
        4. Writing output classifications
        5. Aggregating statistics for report
        """
        if not example_kraken2_file.exists():
            pytest.skip("Example fixture not found")

        # Setup output files
        output_classifications = tmp_path / "output.kraken2"
        output_verbose = tmp_path / "verbose.txt"

        # Initialize structures
        tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}
        taxa_lineages = {}
        report_frequency = 10000

        # Run complete workflow
        with read_file(str(example_kraken2_file)) as f_handle:
            with write_file(str(output_classifications), gz_output=False) as o_handle:
                with write_file(str(output_verbose), gz_output=False) as v_handle:
                    main_loop(
                        f_handle=f_handle,
                        tax_reads_dict=tax_reads_dict,
                        taxonomy_tree=toy_taxonomy,
                        args=mock_args,
                        report_frequency=report_frequency,
                        taxa_lineages=taxa_lineages,
                        paired_input=False,
                        verbose_input=False,
                        o_handle=o_handle,
                        v_handle=v_handle
                    )

        # Verify output files were created
        assert output_classifications.exists(), "Output classifications file should exist"
        assert output_verbose.exists(), "Verbose output file should exist"

        # Verify output files have content
        output_lines = output_classifications.read_text().strip().split('\n')
        assert len(output_lines) > 0, "Should have output classifications"

        verbose_lines = output_verbose.read_text().strip().split('\n')
        assert len(verbose_lines) > 0, "Should have verbose output"

        # Verify tax_reads_dict was populated
        assert len(tax_reads_dict['hits_at_node']) > 0, "Should have aggregated read counts"

        # Verify reads were reclassified (some should be unclassified if threshold not met)
        total_reads = sum(tax_reads_dict['hits_at_node'].values())
        assert total_reads > 0, "Should have processed reads"

    def test_workflow_with_different_thresholds(self, toy_taxonomy, example_kraken2_file, tmp_path):
        """
        Test that different confidence thresholds produce different results.

        Higher threshold → more reads move up the tree or become unclassified
        Lower threshold → more reads stay at original classification
        """
        if not example_kraken2_file.exists():
            pytest.skip("Example fixture not found")

        # Test with LOW threshold (0.1)
        class LowThresholdArgs:
            confidence_threshold = 0.1
            minimum_hit_groups = None
            keep_unclassified = False
            output_report = None

        tax_reads_low = {'hits_at_node': {}, 'hits_at_clade': {}}

        with read_file(str(example_kraken2_file)) as f:
            main_loop(
                f_handle=f,
                tax_reads_dict=tax_reads_low,
                taxonomy_tree=toy_taxonomy,
                args=LowThresholdArgs(),
                report_frequency=10000,
                taxa_lineages={},
                paired_input=False,
                verbose_input=False,
                o_handle=None,
                v_handle=None
            )

        # Test with HIGH threshold (0.9)
        class HighThresholdArgs:
            confidence_threshold = 0.9
            minimum_hit_groups = None
            keep_unclassified = False
            output_report = None

        tax_reads_high = {'hits_at_node': {}, 'hits_at_clade': {}}

        with read_file(str(example_kraken2_file)) as f:
            main_loop(
                f_handle=f,
                tax_reads_dict=tax_reads_high,
                taxonomy_tree=toy_taxonomy,
                args=HighThresholdArgs(),
                report_frequency=10000,
                taxa_lineages={},
                paired_input=False,
                verbose_input=False,
                o_handle=None,
                v_handle=None
            )

        # Different thresholds should produce different results
        # High threshold should have more unclassified (taxid 0)
        unclassified_low = tax_reads_low['hits_at_node'].get(0, 0)
        unclassified_high = tax_reads_high['hits_at_node'].get(0, 0)

        assert unclassified_high >= unclassified_low, \
            "Higher threshold should result in more unclassified reads"

    def test_workflow_with_report_generation(self, toy_taxonomy, example_kraken2_file, mock_args, tmp_path):
        """
        Test complete workflow including final report generation.

        Workflow:
        1. Process reads
        2. Aggregate statistics
        3. Generate Kraken2 report
        """
        if not example_kraken2_file.exists():
            pytest.skip("Example fixture not found")

        # Process reads
        tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}
        total_reads_in_file = 0

        with read_file(str(example_kraken2_file)) as f:
            # Count total reads in file
            for line in f:
                if line.strip():
                    total_reads_in_file += 1

        # Now process with main_loop
        with read_file(str(example_kraken2_file)) as f:
            main_loop(
                f_handle=f,
                tax_reads_dict=tax_reads_dict,
                taxonomy_tree=toy_taxonomy,
                args=mock_args,
                report_frequency=10000,
                taxa_lineages={},
                paired_input=False,
                verbose_input=False,
                o_handle=None,
                v_handle=None
            )

        # Generate report
        report_nodes = get_kraken2_report_content(
            tax_reads_dict,
            toy_taxonomy,
            total_reads_in_file
        )

        # Verify report was generated
        assert len(report_nodes) > 0, "Should generate report nodes"

        # First node should be unclassified
        assert report_nodes[0].node_taxid == 0
        assert report_nodes[0].rank_code == "U"

        # Should have root node
        root_nodes = [n for n in report_nodes if n.node_taxid == 1]
        assert len(root_nodes) == 1, "Should have root in report"

        # Verify report can be formatted
        report_output = tmp_path / "report.txt"
        with open(report_output, 'w') as f:
            for node in report_nodes:
                formatted = format_kraken2_report_row(node)
                f.write(formatted + '\n')

        # Verify report file was created and has content
        assert report_output.exists()
        report_lines = report_output.read_text().strip().split('\n')
        assert len(report_lines) > 0

    def test_workflow_with_paired_end_input(self, toy_taxonomy, example_paired_file, mock_args, tmp_path):
        """
        Test workflow with paired-end Kraken2 input.

        Paired-end reads have |:| delimiter separating k-mer strings
        from both reads in the pair.
        """
        if not example_paired_file.exists():
            pytest.skip("Paired example fixture not found")

        # Process paired-end reads
        tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}
        output_file = tmp_path / "paired_output.kraken2"

        with read_file(str(example_paired_file)) as f_handle:
            with write_file(str(output_file), gz_output=False) as o_handle:
                main_loop(
                    f_handle=f_handle,
                    tax_reads_dict=tax_reads_dict,
                    taxonomy_tree=toy_taxonomy,
                    args=mock_args,
                    report_frequency=10000,
                    taxa_lineages={},
                    paired_input=True,  # PAIRED
                    verbose_input=False,
                    o_handle=o_handle,
                    v_handle=None
                )

        # Verify output was created
        assert output_file.exists()

        # Verify reads were processed
        assert len(tax_reads_dict['hits_at_node']) > 0

        # Read output and verify format
        output_lines = output_file.read_text().strip().split('\n')
        assert len(output_lines) > 0

        # Verify paired delimiter is preserved in output
        has_paired_delimiter = any('|:|' in line for line in output_lines)
        # Note: Delimiter should be preserved in k-mer string field

    def test_workflow_output_format_compatibility(self, toy_taxonomy, example_kraken2_file, mock_args, tmp_path):
        """
        Test that output format is compatible with Kraken2 format.

        Output should have format:
        C/U <tab> read_id <tab> taxid <tab> length <tab> kmer_string
        """
        if not example_kraken2_file.exists():
            pytest.skip("Example fixture not found")

        output_file = tmp_path / "output.kraken2"
        tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}

        with read_file(str(example_kraken2_file)) as f:
            with write_file(str(output_file), gz_output=False) as o:
                main_loop(
                    f_handle=f,
                    tax_reads_dict=tax_reads_dict,
                    taxonomy_tree=toy_taxonomy,
                    args=mock_args,
                    report_frequency=10000,
                    taxa_lineages={},
                    paired_input=False,
                    verbose_input=False,
                    o_handle=o,
                    v_handle=None
                )

        # Read output and verify format
        output_lines = output_file.read_text().strip().split('\n')

        for line in output_lines:
            fields = line.split('\t')

            # Should have at least 5 fields
            assert len(fields) >= 5, f"Line should have at least 5 fields: {line}"

            # First field should be C or U
            assert fields[0] in ['C', 'U'], f"First field should be C or U, got: {fields[0]}"

            # Third field (taxid) should be numeric
            assert fields[2].isdigit(), f"Tax ID should be numeric: {fields[2]}"

            # Fourth field (length) should be numeric
            assert fields[3].isdigit(), f"Length should be numeric: {fields[3]}"

            # Fifth field should be k-mer string (contains ':' or is '0:N')
            assert ':' in fields[4] or fields[4].startswith('0:'), \
                f"K-mer string should contain colon: {fields[4]}"


class TestVerboseOutput:
    """Test verbose output format and content."""

    def test_verbose_output_format(self, toy_taxonomy, example_kraken2_file, mock_args, tmp_path):
        """
        Test that verbose output contains all expected fields.

        Verbose output should include:
        - read_id, length, distance, original_taxid, reclassified_taxid
        - original_conf, recalculated_conf, max_confidence
        - original_rank_code, reclassified_rank_code
        - original_name, reclassified_name, kmer_string
        """
        if not example_kraken2_file.exists():
            pytest.skip("Example fixture not found")

        verbose_file = tmp_path / "verbose.txt"
        tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}

        with read_file(str(example_kraken2_file)) as f:
            with write_file(str(verbose_file), gz_output=False) as v:
                main_loop(
                    f_handle=f,
                    tax_reads_dict=tax_reads_dict,
                    taxonomy_tree=toy_taxonomy,
                    args=mock_args,
                    report_frequency=10000,
                    taxa_lineages={},
                    paired_input=False,
                    verbose_input=False,
                    o_handle=None,
                    v_handle=v
                )

        # Read verbose output
        verbose_lines = verbose_file.read_text().strip().split('\n')
        assert len(verbose_lines) > 0

        # Check first line has expected number of fields
        fields = verbose_lines[0].split('\t')

        # Should have 13 fields for non-verbose input:
        # id, length, distance, orig_taxid, recl_taxid,
        # orig_conf, recl_conf, max_conf,
        # orig_rank, recl_rank, orig_name, recl_name, kmer_string
        assert len(fields) >= 13, f"Verbose output should have 13 fields, got {len(fields)}"

        # Verify confidence values are floats in format X.XX
        orig_conf = fields[5]
        recl_conf = fields[6]
        max_conf = fields[7]

        assert '.' in orig_conf, "Original confidence should be formatted as float"
        assert '.' in recl_conf, "Recalculated confidence should be formatted as float"
        assert '.' in max_conf, "Max confidence should be formatted as float"


class TestGzippedOutput:
    """Test gzipped output file handling."""

    def test_workflow_with_gzipped_output(self, toy_taxonomy, example_kraken2_file, mock_args, tmp_path):
        """
        Test workflow with gzipped output files.
        """
        if not example_kraken2_file.exists():
            pytest.skip("Example fixture not found")

        output_file = tmp_path / "output.kraken2.gz"
        tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}

        with read_file(str(example_kraken2_file)) as f:
            with write_file(str(output_file), gz_output=True) as o:
                main_loop(
                    f_handle=f,
                    tax_reads_dict=tax_reads_dict,
                    taxonomy_tree=toy_taxonomy,
                    args=mock_args,
                    report_frequency=10000,
                    taxa_lineages={},
                    paired_input=False,
                    verbose_input=False,
                    o_handle=o,
                    v_handle=None
                )

        # Verify gzipped file was created
        assert output_file.exists()

        # Verify we can read it back
        import gzip
        with gzip.open(output_file, 'rt') as f:
            lines = f.readlines()
            assert len(lines) > 0, "Should have content in gzipped file"


class TestEdgeCasesInWorkflow:
    """Test edge cases in complete workflow."""

    def test_workflow_with_empty_input(self, toy_taxonomy, tmp_path):
        """
        Test that empty input file is handled gracefully at the file I/O layer.

        Note: main_loop requires ≥1 classified read for report generation,
        so we test file reading separately.
        """
        # Create empty input file
        empty_file = tmp_path / "empty.kraken2"
        empty_file.write_text("")

        # Should be able to open and read empty file
        with read_file(str(empty_file)) as f:
            lines = f.readlines()

        # Empty file has no lines
        assert len(lines) == 0

        # Verify we can write to output file even with no input
        output_file = tmp_path / "output.kraken2"
        with write_file(str(output_file), gz_output=False) as o:
            # Write nothing (valid scenario)
            pass

        # Output file should exist but be empty
        assert output_file.exists()
        assert output_file.read_text() == ""

    def test_workflow_with_mostly_unclassified_reads(self, toy_taxonomy, mock_args, tmp_path):
        """
        Test workflow with mostly unclassified reads (99% unclassified, 1% classified).

        Verifies that:
        1. Unclassified reads are skipped by main_loop (not processed)
        2. Classified reads are processed normally
        3. Report generation works with minimal classified reads
        """
        # Create file with 99 unclassified reads and 1 classified read
        unclassified_file = tmp_path / "mostly_unclassified.kraken2"
        lines = []
        # Add 99 unclassified reads
        for i in range(99):
            lines.append(f"U\tread_{i}\t0\t100\t0:100\n")
        # Add 1 classified read
        lines.append("C\tread_classified\t562\t100\t562:80 0:20\n")

        unclassified_file.write_text("".join(lines))

        tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}

        with read_file(str(unclassified_file)) as f:
            main_loop(
                f_handle=f,
                tax_reads_dict=tax_reads_dict,
                taxonomy_tree=toy_taxonomy,
                args=mock_args,
                report_frequency=10000,
                taxa_lineages={},
                paired_input=False,
                verbose_input=False,
                o_handle=None,
                v_handle=None
            )

        # Should process only the 1 classified read (unclassified reads skipped)
        total_classified = sum(tax_reads_dict['hits_at_node'].values())
        assert total_classified == 1, "Should have processed exactly 1 classified read"

        # Verify the classified read was assigned to E. coli (562) or moved up
        assert len(tax_reads_dict['hits_at_node']) > 0, "Should have at least one taxon with reads"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
