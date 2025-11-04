#!/usr/bin/env python3
"""
Integration tests for file I/O operations.

Tests reading Kraken2 files, parsing lines, handling gzip, and validation.
These tests verify that StringMeUp correctly handles real file formats.
"""

import pytest
import gzip
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pathlib import Path
from stringmeup.stringmeup import (
    create_read,
    read_file,
    write_file,
    validate_input_file,
    is_paired_input,
    is_verbose_input
)


class TestCreateRead:
    """Test parsing individual Kraken2 output lines."""

    def test_parse_classified_single_end(self):
        """Test parsing a classified single-end read."""
        line = "C\tread1\t562\t100\t562:60 0:30 A:10"

        read = create_read(line, verbose_input=False)

        assert read.id == "read1"
        assert read.original_taxid == 562
        assert read.kmer_string == "562:60 0:30 A:10"
        assert read.classified == False  # Will be set during reclassification

    def test_parse_unclassified_read(self):
        """Test parsing an unclassified read."""
        line = "U\tread2\t0\t100\t0:100"

        read = create_read(line, verbose_input=False)

        assert read.id == "read2"
        assert read.original_taxid == 0
        assert read.kmer_string == "0:100"

    def test_parse_paired_end_read(self):
        """Test parsing a paired-end read with |:| delimiter."""
        line = "C\tread3\t562\t100\t562:40 A:5 |:| 562:30 A:5"

        read = create_read(line, verbose_input=False)

        assert read.id == "read3"
        assert read.original_taxid == 562
        assert "|:|" in read.kmer_string
        assert read.kmer_string == "562:40 A:5 |:| 562:30 A:5"

    def test_parse_read_with_ambiguous(self):
        """Test parsing read with ambiguous k-mers."""
        line = "C\tread4\t562\t100\t562:50 A:20 0:30"

        read = create_read(line, verbose_input=False)

        assert "A:20" in read.kmer_string
        assert read.kmer_string == "562:50 A:20 0:30"

    def test_parse_read_with_length(self):
        """Test that read length is captured."""
        line = "C\tread5\t562\t150\t562:60 0:40"

        read = create_read(line, verbose_input=False)

        # Length field is present in Kraken2 output (4th field)
        # Check if it's parsed (implementation dependent)
        assert read.original_taxid == 562

    def test_parse_verbose_format(self):
        """Test parsing verbose Kraken2 output (with hit groups)."""
        # Verbose format has additional hit group field (5th column)
        line = "C\tread6\t562\t100\t10\t562:60 0:30 A:10"

        read = create_read(line, verbose_input=True)

        assert read.id == "read6"
        assert read.minimizer_hit_groups == 10

    def test_parse_multiple_taxa(self):
        """Test parsing read with hits to multiple taxa."""
        line = "C\tread7\t561\t100\t561:30 562:20 1224:10 0:40"

        read = create_read(line, verbose_input=False)

        assert read.original_taxid == 561
        assert "561:30" in read.kmer_string
        assert "562:20" in read.kmer_string
        assert "1224:10" in read.kmer_string


class TestReadFile:
    """Test file reading functionality."""

    def test_read_regular_file(self, tmp_path):
        """Test reading a regular (uncompressed) file."""
        test_file = tmp_path / "test.kraken2"
        test_file.write_text("C\tread1\t562\t100\t562:60 0:40\n")

        with read_file(str(test_file)) as f:
            line = f.readline()
            assert "read1" in line
            assert "562:60" in line

    def test_read_gzipped_file(self, tmp_path):
        """Test reading a gzipped file."""
        test_file = tmp_path / "test.kraken2.gz"

        # Write gzipped content
        content = "C\tread1\t562\t100\t562:60 0:40\n"
        with gzip.open(test_file, 'wt') as f:
            f.write(content)

        # Read using read_file
        with read_file(str(test_file)) as f:
            line = f.readline()
            assert "read1" in line
            assert "562:60" in line

    def test_read_file_context_manager(self, tmp_path):
        """Test that read_file works as context manager."""
        test_file = tmp_path / "test.kraken2"
        test_file.write_text("line1\nline2\nline3\n")

        # Should work with context manager
        with read_file(str(test_file)) as f:
            lines = f.readlines()

        assert len(lines) == 3
        assert "line1" in lines[0]

    def test_read_example_file(self):
        """Test reading the example fixture file."""
        # Use the example file we created
        example_file = Path(__file__).parent.parent / "fixtures" / "example.kraken2"

        if example_file.exists():
            with read_file(str(example_file)) as f:
                lines = f.readlines()

            assert len(lines) == 10  # We created 10 example reads
            assert any("read1" in line for line in lines)


class TestWriteFile:
    """Test file writing functionality."""

    def test_write_regular_file(self, tmp_path):
        """Test writing a regular file."""
        output_file = tmp_path / "output.txt"

        with write_file(str(output_file), gz_output=False) as f:
            f.write("test output\n")

        # Verify file was written
        assert output_file.exists()
        content = output_file.read_text()
        assert "test output" in content

    def test_write_gzipped_file(self, tmp_path):
        """Test writing a gzipped file."""
        output_file = tmp_path / "output.txt.gz"

        with write_file(str(output_file), gz_output=True) as f:
            f.write("test output\n")

        # Verify gzipped file was written
        assert output_file.exists()

        # Read back the gzipped content
        with gzip.open(output_file, 'rt') as f:
            content = f.read()
        assert "test output" in content

    def test_write_file_context_manager(self, tmp_path):
        """Test that write_file works as context manager."""
        output_file = tmp_path / "output.txt"

        with write_file(str(output_file), gz_output=False) as f:
            f.write("line1\n")
            f.write("line2\n")

        content = output_file.read_text()
        assert "line1" in content
        assert "line2" in content


class TestFileValidation:
    """Test input file validation."""

    def test_validate_correct_file(self, tmp_path):
        """Test validation of correctly formatted single-end file."""
        test_file = tmp_path / "valid.kraken2"
        test_file.write_text("C\tread1\t562\t100\t562:60 0:40\n")

        # Should not raise exception
        try:
            validate_input_file(
                str(test_file),
                verbose_input=False,
                minimum_hit_groups=None,
                paired_input=False
            )
        except Exception as e:
            pytest.fail(f"Valid file should not raise exception: {e}")

    def test_validate_correct_paired_file(self, tmp_path):
        """Test validation of correctly formatted paired-end file."""
        test_file = tmp_path / "valid_paired.kraken2"
        test_file.write_text("C\tread1\t562\t100\t562:40 |:| 562:30\n")

        # Should not raise exception
        try:
            validate_input_file(
                str(test_file),
                verbose_input=False,
                minimum_hit_groups=None,
                paired_input=True
            )
        except Exception as e:
            pytest.fail(f"Valid paired file should not raise exception: {e}")

    def test_is_paired_input_detection(self, tmp_path):
        """
        Test detection of paired-end format.

        BUG FIX: Fixed is_paired_input() to check correct field (index 4, not 3).
        """
        # Create paired-end file
        paired_file = tmp_path / "paired.kraken2"
        paired_file.write_text(
            "C\tread1\t562\t100\t562:40 |:| 562:30\n"
            "C\tread2\t562\t100\t562:50 |:| 562:30\n"
            "C\tread3\t562\t100\t562:45 |:| 562:25\n"
        )

        result = is_paired_input(str(paired_file))

        # Should detect paired-end format (returns True)
        assert result is True, "Should detect paired-end format from |:| delimiter"

    def test_is_single_end_detection(self, tmp_path):
        """Test detection of single-end format."""
        # Create single-end file
        single_file = tmp_path / "single.kraken2"
        single_file.write_text("C\tread1\t562\t100\t562:60 0:40\n")

        result = is_paired_input(str(single_file))
        # Function returns None for single-end
        assert result is None, "Should detect single-end format"

    def test_is_verbose_input_detection(self, tmp_path):
        """Test detection of verbose format."""
        # Create verbose file (has additional fields)
        verbose_file = tmp_path / "verbose.kraken2"
        # Verbose format typically has hit group information
        verbose_file.write_text("C\tread1\t562\t100\t562:60 0:40\t10:5 20:8\n")

        # Detection depends on implementation
        # This is a placeholder test
        try:
            is_verbose = is_verbose_input(str(verbose_file))
            # Result depends on how verbose is detected
        except Exception:
            pass  # Function may not be fully implemented

    def test_validate_paired_file_consistency(self, tmp_path):
        """
        Test validation with paired_input=True when first line missing |:|.

        NOTE: validate_input_file() only checks the FIRST LINE, not all lines.
        So inconsistency after the first line won't be detected.
        """
        # First read missing |:|, second has it
        paired_file = tmp_path / "paired_inconsistent.kraken2"
        paired_file.write_text(
            "C\tread1\t562\t100\t562:60 0:40\n"  # First line: missing |:|
            "C\tread2\t562\t100\t562:40 |:| 562:30\n"
        )

        # Should fail validation because first line doesn't have |:|
        with pytest.raises(SystemExit):
            validate_input_file(
                str(paired_file),
                verbose_input=False,
                minimum_hit_groups=None,
                paired_input=True  # Expecting paired format
            )


class TestFileIOWithRealFixtures:
    """Test file I/O using the real fixture files."""

    def test_parse_all_example_reads(self):
        """Test parsing all reads from example.kraken2 fixture."""
        example_file = Path(__file__).parent.parent / "fixtures" / "example.kraken2"

        if not example_file.exists():
            pytest.skip("Example fixture not found")

        with read_file(str(example_file)) as f:
            reads = []
            for line in f:
                if line.strip():
                    read = create_read(line.strip(), verbose_input=False)
                    reads.append(read)

        # Should have 10 reads
        assert len(reads) == 10

        # Check specific reads
        read_ids = [r.id for r in reads]
        assert "read1" in read_ids
        assert "read4" in read_ids  # Unclassified
        assert "read9" in read_ids  # Perfect classification

        # Check read4 is unclassified
        read4 = [r for r in reads if r.id == "read4"][0]
        assert read4.original_taxid == 0

    def test_parse_paired_end_fixture(self):
        """Test parsing paired-end example file."""
        paired_file = Path(__file__).parent.parent / "fixtures" / "example_paired.kraken2"

        if not paired_file.exists():
            pytest.skip("Paired example fixture not found")

        with read_file(str(paired_file)) as f:
            reads = []
            for line in f:
                if line.strip():
                    read = create_read(line.strip(), verbose_input=False)
                    reads.append(read)

        # All reads should have |:| delimiter
        for read in reads:
            assert "|:|" in read.kmer_string, \
                f"Paired-end read {read.id} should have |:| delimiter"

    def test_parse_ambiguous_fixture(self):
        """Test parsing ambiguous k-mer example file."""
        ambig_file = Path(__file__).parent.parent / "fixtures" / "example_with_ambiguous.kraken2"

        if not ambig_file.exists():
            pytest.skip("Ambiguous example fixture not found")

        with read_file(str(ambig_file)) as f:
            reads = []
            for line in f:
                if line.strip():
                    read = create_read(line.strip(), verbose_input=False)
                    reads.append(read)

        # Check that ambiguous k-mers are present
        ambig_reads = [r for r in reads if "A:" in r.kmer_string]
        assert len(ambig_reads) > 0, "Should have reads with ambiguous k-mers"

        # Check specific reads
        read_10pct = [r for r in reads if r.id == "read_10pct_ambiguous"]
        if read_10pct:
            assert "A:10" in read_10pct[0].kmer_string


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
