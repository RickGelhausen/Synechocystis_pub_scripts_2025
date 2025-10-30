"""
Unit tests for lib.bam_reader module

Tests for:
- LocusExtractor class
- IntervalReader class
"""

import pytest
import sys
from pathlib import Path
from unittest.mock import Mock, patch
import interlap


from lib.bam_reader import LocusExtractor, IntervalReader


class TestLocusExtractor:
    """Test suite for LocusExtractor class"""

    @pytest.fixture
    def mock_alignment_file(self):
        """Create a mock alignment file path"""
        mock_path = Mock(spec=Path)
        mock_path.stem = "test_sample"
        mock_path.__fspath__ = Mock(return_value="/path/to/test_sample.bam")
        mock_path.__str__ = Mock(return_value="/path/to/test_sample.bam")
        return mock_path

    @pytest.fixture
    def config(self):
        """Standard configuration for tests"""
        return {
            "positionsOutsideORF": 50,
            "positionsInsideORF": 100
        }

    @pytest.fixture
    def empty_dicts(self):
        """Empty output dictionaries"""
        return {}, {}

    @pytest.fixture
    def start_loci_interlaps(self):
        """Mock start loci interlaps"""
        interlaps = {}
        inter = interlap.InterLap()
        # Add some test intervals: (start, stop, gene_id)
        inter.update([(100, 200, "gene1"), (300, 400, "gene2")])
        interlaps[("chr1", "+")] = inter

        inter_minus = interlap.InterLap()
        inter_minus.update([(500, 600, "gene3")])
        interlaps[("chr1", "-")] = inter_minus

        return interlaps

    @pytest.fixture
    def stop_loci_interlaps(self):
        """Mock stop loci interlaps"""
        interlaps = {}
        inter = interlap.InterLap()
        inter.update([(150, 250, "gene1"), (350, 450, "gene2")])
        interlaps[("chr1", "+")] = inter
        return interlaps

    @pytest.fixture
    def mock_reads(self):
        """Create mock reads for testing"""
        reads = []

        # Read 1: Plus strand, length 30
        read1 = Mock()
        read1.reference_name = "chr1"
        read1.get_tag.return_value = 1  # NH tag
        read1.mapping_quality = 42
        read1.is_unmapped = False
        read1.reference_start = 150
        read1.query_length = 30
        read1.is_reverse = False
        reads.append(read1)

        # Read 2: Minus strand, length 32
        read2 = Mock()
        read2.reference_name = "chr1"
        read2.get_tag.return_value = 1
        read2.mapping_quality = 42
        read2.is_unmapped = False
        read2.reference_start = 550
        read2.query_length = 32
        read2.is_reverse = True
        reads.append(read2)

        # Read 3: Should be filtered (NH > 1)
        read3 = Mock()
        read3.reference_name = "chr1"
        read3.get_tag.return_value = 2  # Multiple alignments
        read3.mapping_quality = 42
        read3.is_unmapped = False
        reads.append(read3)

        # Read 4: Should be filtered (unmapped)
        read4 = Mock()
        read4.reference_name = "chr1"
        read4.get_tag.return_value = 1
        read4.mapping_quality = 42
        read4.is_unmapped = True
        reads.append(read4)

        return reads

    @patch('pysam.AlignmentFile')
    def test_init_calls_extract_loci(self, mock_pysam, mock_alignment_file,
                                     start_loci_interlaps, stop_loci_interlaps,
                                     config, empty_dicts):
        """Test that initialization calls _extract_loci"""
        mock_bam = Mock()
        mock_bam.fetch.return_value = []
        mock_pysam.return_value = mock_bam

        start_dict, stop_dict = empty_dicts

        with patch('builtins.print'):
            _ = LocusExtractor(
                mock_alignment_file,
                "threeprime",
                start_loci_interlaps,
                stop_loci_interlaps,
                config,
                start_dict,
                stop_dict
            )

        mock_pysam.assert_called_once_with(mock_alignment_file)
        mock_bam.fetch.assert_called_once()

    def test_init_output_dict_creates_nested_structure(self, mock_alignment_file,
                                                        start_loci_interlaps,
                                                        stop_loci_interlaps,
                                                        config, empty_dicts):
        """Test _init_output_dict creates proper nested dictionary structure"""
        start_dict, stop_dict = empty_dicts

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_bam = Mock()
            mock_bam.fetch.return_value = []
            mock_pysam.return_value = mock_bam

            with patch('builtins.print'):
                extractor = LocusExtractor(
                    mock_alignment_file,
                    "threeprime",
                    start_loci_interlaps,
                    stop_loci_interlaps,
                    config,
                    start_dict,
                    stop_dict
                )

            output_dict = {}
            extractor._init_output_dict(output_dict, "chr1", 30)

            assert "chr1" in output_dict
            assert "test_sample" in output_dict["chr1"]
            assert 30 in output_dict["chr1"]["test_sample"]
            assert isinstance(output_dict["chr1"]["test_sample"][30], dict)

    @pytest.mark.parametrize("strand,mapping_mode,start,stop,expected", [
        ("+", "threeprime", 100, 129, 129),  # Plus strand, 3' end
        ("+", "fiveprime", 100, 129, 100),   # Plus strand, 5' end
        ("-", "threeprime", 100, 129, 100),  # Minus strand, 3' end
        ("-", "fiveprime", 100, 129, 129),   # Minus strand, 5' end
    ])
    def test_get_genome_position(self, mock_alignment_file, start_loci_interlaps,
                                 stop_loci_interlaps, config, empty_dicts,
                                 strand, mapping_mode, start, stop, expected):
        """Test genome position calculation for different strands and modes"""
        start_dict, stop_dict = empty_dicts

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_bam = Mock()
            mock_bam.fetch.return_value = []
            mock_pysam.return_value = mock_bam

            with patch('builtins.print'):
                extractor = LocusExtractor(
                    mock_alignment_file,
                    mapping_mode,
                    start_loci_interlaps,
                    stop_loci_interlaps,
                    config,
                    start_dict,
                    stop_dict
                )

            result = extractor._get_genome_position(strand, start, stop)
            assert result == expected

    def test_get_genome_position_invalid_mode(self, mock_alignment_file,
                                               start_loci_interlaps,
                                               stop_loci_interlaps,
                                               config, empty_dicts):
        """Test that invalid mapping mode raises SystemExit"""
        start_dict, stop_dict = empty_dicts

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_bam = Mock()
            mock_bam.fetch.return_value = []
            mock_pysam.return_value = mock_bam

            with patch('builtins.print'):
                extractor = LocusExtractor(
                    mock_alignment_file,
                    "invalid_mode",
                    start_loci_interlaps,
                    stop_loci_interlaps,
                    config,
                    start_dict,
                    stop_dict
                )

            with pytest.raises(SystemExit):
                extractor._get_genome_position("+", 100, 129)

    @pytest.mark.parametrize("strand,position_genome,window,interval,expected", [
        ("+", 150, [-50, 100], (100, 250, "start_window1"), 0),   # Plus strand at interval start
        ("+", 175, [-50, 100], (100, 250, "start_window1"), 25),  # Plus strand middle
        ("-", 150, [-50, 100], (100, 250, "start_window1"), 50),  # Minus strand at interval start
        ("-", 175, [-50, 100], (100, 250, "start_window1"), 25),  # Minus strand middle
    ])
    def test_calculate_window_position(self, mock_alignment_file,
                                       start_loci_interlaps,
                                       stop_loci_interlaps,
                                       config, empty_dicts,
                                       strand, position_genome, window,
                                       interval, expected):
        """Test window position calculation relative to interval"""
        start_dict, stop_dict = empty_dicts

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_bam = Mock()
            mock_bam.fetch.return_value = []
            mock_pysam.return_value = mock_bam

            with patch('builtins.print'):
                extractor = LocusExtractor(
                    mock_alignment_file,
                    "threeprime",
                    start_loci_interlaps,
                    stop_loci_interlaps,
                    config,
                    start_dict,
                    stop_dict
                )

            result = extractor._calculate_window_position(
                strand, position_genome, window, interval
            )
            assert result == expected

    @patch('pysam.AlignmentFile')
    def test_extract_loci_filters_multimapped_reads(self, mock_pysam,
                                                     mock_alignment_file,
                                                     start_loci_interlaps,
                                                     stop_loci_interlaps,
                                                     config, empty_dicts,
                                                     mock_reads):
        """Test that reads with NH > 1 are filtered out"""
        mock_bam = Mock()
        mock_bam.fetch.return_value = mock_reads
        mock_pysam.return_value = mock_bam

        start_dict, stop_dict = empty_dicts

        with patch('builtins.print'):
            _ = LocusExtractor(
                mock_alignment_file,
                "threeprime",
                start_loci_interlaps,
                stop_loci_interlaps,
                config,
                start_dict,
                stop_dict
            )

        # Only 2 reads should be processed (read1 and read2)
        # read3 (NH=2) and read4 (unmapped) should be filtered
        assert len(start_dict.get("chr1", {}).get("test_sample", {})) == 2

    @patch('pysam.AlignmentFile')
    def test_extract_loci_processes_valid_reads(self, mock_pysam,
                                                mock_alignment_file,
                                                start_loci_interlaps,
                                                stop_loci_interlaps,
                                                config, empty_dicts,
                                                mock_reads):
        """Test that valid reads are properly processed and counted"""
        mock_bam = Mock()
        mock_bam.fetch.return_value = [mock_reads[0]]  # Just the first valid read
        mock_pysam.return_value = mock_bam

        start_dict, stop_dict = empty_dicts

        with patch('builtins.print'):
            extractor = LocusExtractor(
                mock_alignment_file,
                "threeprime",
                start_loci_interlaps,
                stop_loci_interlaps,
                config,
                start_dict,
                stop_dict
            )

        # Check that dictionaries were populated
        assert "chr1" in start_dict
        assert "test_sample" in start_dict["chr1"]

    @patch('pysam.AlignmentFile')
    def test_extract_loci_handles_missing_index(self, mock_pysam,
                                                mock_alignment_file,
                                                start_loci_interlaps,
                                                stop_loci_interlaps,
                                                config, empty_dicts):
        """Test that ValueError from missing index is handled"""
        mock_bam = Mock()
        mock_bam.fetch.side_effect = ValueError("Index not found")
        mock_pysam.return_value = mock_bam

        start_dict, stop_dict = empty_dicts

        with patch('builtins.print'):
            with pytest.raises(SystemExit):
                extractor = LocusExtractor(
                    mock_alignment_file,
                    "threeprime",
                    start_loci_interlaps,
                    stop_loci_interlaps,
                    config,
                    start_dict,
                    stop_dict
                )

    @patch('pysam.AlignmentFile')
    def test_output_returns_dictionaries(self, mock_pysam, mock_alignment_file,
                                        start_loci_interlaps, stop_loci_interlaps,
                                        config, empty_dicts):
        """Test that output() returns both dictionaries"""
        mock_bam = Mock()
        mock_bam.fetch.return_value = []
        mock_pysam.return_value = mock_bam

        start_dict, stop_dict = empty_dicts

        with patch('builtins.print'):
            extractor = LocusExtractor(
                mock_alignment_file,
                "threeprime",
                start_loci_interlaps,
                stop_loci_interlaps,
                config,
                start_dict,
                stop_dict
            )

        result_start, result_stop = extractor.output()
        assert result_start is start_dict
        assert result_stop is stop_dict


# class TestIntervalReader:
#     """Test suite for IntervalReader class"""

#     @pytest.fixture
#     def mock_alignment_file(self):
#         """Create a mock alignment file path"""
#         mock_path = Mock(spec=Path)
#         mock_path.stem = "test_intervals"
#         mock_path.__str__ = lambda x: "/path/to/test_intervals.bam"
#         return mock_path

#     @pytest.fixture
#     def mock_reads(self):
#         """Create mock reads for IntervalReader testing"""
#         reads = []

#         # Read 1: chr1, plus strand
#         read1 = Mock()
#         read1.reference_name = "chr1"
#         read1.get_tag.return_value = 1
#         read1.mapping_quality = 42
#         read1.is_unmapped = False
#         read1.reference_start = 100
#         read1.query_length = 30
#         read1.is_reverse = False
#         reads.append(read1)

#         # Read 2: chr1, plus strand (another read on same chromosome)
#         read2 = Mock()
#         read2.reference_name = "chr1"
#         read2.get_tag.return_value = 1
#         read2.mapping_quality = 42
#         read2.is_unmapped = False
#         read2.reference_start = 200
#         read2.query_length = 32
#         read2.is_reverse = False
#         reads.append(read2)

#         # Read 3: chr1, minus strand
#         read3 = Mock()
#         read3.reference_name = "chr1"
#         read3.get_tag.return_value = 1
#         read3.mapping_quality = 42
#         read3.is_unmapped = False
#         read3.reference_start = 300
#         read3.query_length = 28
#         read3.is_reverse = True
#         reads.append(read3)

#         # Read 4: chr2, plus strand
#         read4 = Mock()
#         read4.reference_name = "chr2"
#         read4.get_tag.return_value = 1
#         read4.mapping_quality = 42
#         read4.is_unmapped = False
#         read4.reference_start = 150
#         read4.query_length = 30
#         read4.is_reverse = False
#         reads.append(read4)

#         # Read 5: Should be filtered (multimapped)
#         read5 = Mock()
#         read5.reference_name = "chr1"
#         read5.get_tag.return_value = 3
#         read5.mapping_quality = 0
#         read5.is_unmapped = False
#         reads.append(read5)

#         return reads

#     @patch('pysam.AlignmentFile')
#     def test_init_calls_read_alignment_file(self, mock_pysam,
#                                            mock_alignment_file, mock_reads):
#         """Test that initialization calls _read_alignment_file"""
#         mock_bam = Mock()
#         mock_bam.fetch.return_value = mock_reads
#         mock_pysam.return_value = mock_bam

#         with patch('builtins.print'):
#             reader = IntervalReader(mock_alignment_file)

#         mock_pysam.assert_called_once_with(mock_alignment_file)
#         mock_bam.fetch.assert_called_once()

#     @patch('pysam.AlignmentFile')
#     def test_read_alignment_file_counts_reads_per_chromosome(self, mock_pysam,
#                                                              mock_alignment_file,
#                                                              mock_reads):
#         """Test that reads are properly counted per chromosome"""
#         mock_bam = Mock()
#         mock_bam.fetch.return_value = mock_reads[:4]  # Exclude filtered read
#         mock_pysam.return_value = mock_bam

#         with patch('builtins.print'):
#             reader = IntervalReader(mock_alignment_file)

#         # 3 reads on chr1, 1 read on chr2
#         assert reader.no_accepted_reads_dict["chr1"] == 3
#         assert reader.no_accepted_reads_dict["chr2"] == 1

#     @patch('pysam.AlignmentFile')
#     def test_read_alignment_file_creates_interlaps_by_strand(self, mock_pysam,
#                                                              mock_alignment_file,
#                                                              mock_reads):
#         """Test that InterLap objects are created for each chromosome/strand combination"""
#         mock_bam = Mock()
#         mock_bam.fetch.return_value = mock_reads[:4]
#         mock_pysam.return_value = mock_bam

#         with patch('builtins.print'):
#             reader = IntervalReader(mock_alignment_file)

#         # Should have 3 combinations: (chr1, +), (chr1, -), (chr2, +)
#         assert ("chr1", "+") in reader.reads_interlap_dict
#         assert ("chr1", "-") in reader.reads_interlap_dict
#         assert ("chr2", "+") in reader.reads_interlap_dict
#         assert ("chr2", "-") not in reader.reads_interlap_dict

#     @patch('pysam.AlignmentFile')
#     def test_read_alignment_file_filters_invalid_reads(self, mock_pysam,
#                                                        mock_alignment_file,
#                                                        mock_reads):
#         """Test that reads with NH > 1, low quality, or unmapped are filtered"""
#         mock_bam = Mock()
#         mock_bam.fetch.return_value = mock_reads  # Includes filtered reads
#         mock_pysam.return_value = mock_bam

#         with patch('builtins.print'):
#             reader = IntervalReader(mock_alignment_file)

#         # Should only count the 4 valid reads, not the 5th one
#         total_reads = sum(reader.no_accepted_reads_dict.values())
#         assert total_reads == 4

#     @patch('pysam.AlignmentFile')
#     def test_read_alignment_file_stores_read_length_in_intervals(self, mock_pysam,
#                                                                  mock_alignment_file,
#                                                                  mock_reads):
#         """Test that intervals store read length information"""
#         mock_bam = Mock()
#         mock_bam.fetch.return_value = [mock_reads[0]]  # Single read
#         mock_pysam.return_value = mock_bam

#         with patch('builtins.print'):
#             reader = IntervalReader(mock_alignment_file)

#         interlap_obj = reader.reads_interlap_dict[("chr1", "+")]
#         # Find intervals overlapping the read position
#         intervals = list(interlap_obj.find((100, 129)))

#         assert len(intervals) == 1
#         start, stop, read_length = intervals[0]
#         assert start == 100
#         assert stop == 129  # 100 + 30 - 1
#         assert read_length == 30

#     @patch('pysam.AlignmentFile')
#     def test_read_alignment_file_handles_missing_index(self, mock_pysam,
#                                                        mock_alignment_file):
#         """Test that ValueError from missing index is handled"""
#         mock_bam = Mock()
#         mock_bam.fetch.side_effect = ValueError("Index not found")
#         mock_pysam.return_value = mock_bam

#         with patch('builtins.print'):
#             with pytest.raises(SystemExit):
#                 reader = IntervalReader(mock_alignment_file)

#     @patch('pysam.AlignmentFile')
#     def test_output_returns_correct_data_structures(self, mock_pysam,
#                                                     mock_alignment_file,
#                                                     mock_reads):
#         """Test that output() returns both dictionaries with correct structure"""
#         mock_bam = Mock()
#         mock_bam.fetch.return_value = mock_reads[:4]
#         mock_pysam.return_value = mock_bam

#         with patch('builtins.print'):
#             reader = IntervalReader(mock_alignment_file)

#         interlaps, read_counts = reader.output()

#         # Check interlaps structure
#         assert isinstance(interlaps, dict)
#         for key, value in interlaps.items():
#             assert isinstance(key, tuple)
#             assert len(key) == 2  # (chrom, strand)
#             assert isinstance(value, interlap.InterLap)

#         # Check read counts structure
#         assert isinstance(read_counts, dict)
#         assert "chr1" in read_counts
#         assert "chr2" in read_counts

#     @patch('pysam.AlignmentFile')
#     def test_empty_bam_file(self, mock_pysam, mock_alignment_file):
#         """Test handling of empty BAM file"""
#         mock_bam = Mock()
#         mock_bam.fetch.return_value = []
#         mock_pysam.return_value = mock_bam

#         with patch('builtins.print'):
#             reader = IntervalReader(mock_alignment_file)

#         interlaps, read_counts = reader.output()

#         assert len(interlaps) == 0
#         assert len(read_counts) == 0


# # Integration-style tests
# class TestLocusExtractorIntegration:
#     """Integration tests for LocusExtractor with more realistic scenarios"""

#     @patch('pysam.AlignmentFile')
#     def test_full_workflow_threeprime_plus_strand(self, mock_pysam):
#         """Test complete workflow for 3' mapping on plus strand"""
#         # Setup
#         mock_path = Mock(spec=Path)
#         mock_path.stem = "sample1"

#         # Create a read that overlaps with a gene
#         read = Mock()
#         read.reference_name = "chr1"
#         read.get_tag.return_value = 1
#         read.mapping_quality = 42
#         read.is_unmapped = False
#         read.reference_start = 140
#         read.query_length = 30
#         read.is_reverse = False

#         mock_bam = Mock()
#         mock_bam.fetch.return_value = [read]
#         mock_pysam.return_value = mock_bam

#         # Create interlaps with a gene
#         start_interlaps = {}
#         inter = interlap.InterLap()
#         inter.update([(100, 200, "gene1")])
#         start_interlaps[("chr1", "+")] = inter

#         stop_interlaps = {}

#         config = {
#             "positionsOutsideORF": 50,
#             "positionsInsideORF": 100
#         }

#         start_dict = {}
#         stop_dict = {}

#         with patch('builtins.print'):
#             extractor = LocusExtractor(
#                 mock_path,
#                 "threeprime",
#                 start_interlaps,
#                 stop_interlaps,
#                 config,
#                 start_dict,
#                 stop_dict
#             )

#         # Verify results
#         assert "chr1" in start_dict
#         assert "sample1" in start_dict["chr1"]
#         assert 30 in start_dict["chr1"]["sample1"]

#         # Check that position was mapped correctly
#         positions = list(start_dict["chr1"]["sample1"][30].keys())
#         assert len(positions) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])