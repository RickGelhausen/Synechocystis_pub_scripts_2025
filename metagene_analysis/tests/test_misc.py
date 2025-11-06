
from lib.misc import (
    calculate_rpkm,
    count_reads,
    window_normalize_df,
    json_dict_to_dataframe,
    aggregate_coverage_dataframe,
    apply_normalization
)

import pytest
import interlap
import pandas as pd
import numpy as np

class TestCalculateRpkm:
    """Test suite for calculate_rpkm function"""

    def test_basic_calculation(self):
        """Test basic RPKM calculation"""
        gene_length = 1000
        read_counts = 100
        total_counts = 10000000

        result = calculate_rpkm(gene_length, read_counts, total_counts)

        # RPKM = (reads * 1e9) / (total_counts * gene_length)
        # RPKM = (100 * 1e9) / (10000000 * 1000)
        # RPKM = 1e11 / 1e10 = 10.0
        assert result == 10.0

    def test_zero_reads(self):
        """Test with zero read counts"""
        result = calculate_rpkm(1000, 0, 10000000)
        assert result == 0.0

    def test_high_expression(self):
        """Test with high expression levels"""
        gene_length = 500
        read_counts = 5000
        total_counts = 1000000

        result = calculate_rpkm(gene_length, read_counts, total_counts)
        expected = (5000 * 1e9) / (1000000 * 500)

        assert result == expected

    def test_low_expression(self):
        """Test with low expression levels"""
        gene_length = 2000
        read_counts = 1
        total_counts = 100000000

        result = calculate_rpkm(gene_length, read_counts, total_counts)
        expected = (1 * 1e9) / (100000000 * 2000)

        assert result == expected

    def test_returns_float(self):
        """Test that result is a float"""
        result = calculate_rpkm(1000, 100, 10000000)
        assert isinstance(result, float)



class TestCountReads:
    """Test suite for count_reads function"""

    @pytest.fixture
    def mock_read_intervals_dict(self):
        """Create mock read intervals dictionary"""
        intervals_dict = {}

        # Plus strand intervals
        inter_plus = interlap.InterLap()
        inter_plus.update([
            (100, 129, 30),  # Read length 30
            (150, 179, 30),
            (200, 229, 30),
        ])
        intervals_dict[("chr1", "+")] = inter_plus

        # Minus strand intervals
        inter_minus = interlap.InterLap()
        inter_minus.update([
            (100, 129, 30),
            (150, 179, 30),
        ])
        intervals_dict[("chr1", "-")] = inter_minus

        return intervals_dict

    def test_global_mapping_counts_all(self, mock_read_intervals_dict):
        """Test global mapping counts all overlapping reads"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 100, 250, "global"
        )
        # All 3 reads overlap the region
        assert count == 3

    def test_threeprime_plus_strand(self, mock_read_intervals_dict):
        """Test 3' mapping on plus strand (uses read end)"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 100, 180, "threeprime"
        )
        # Reads with end (position [1]) <= 180: (100,129) and (150,179)
        assert count == 2

    def test_threeprime_minus_strand(self, mock_read_intervals_dict):
        """Test 3' mapping on minus strand (uses read start)"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "-", 100, 180, "threeprime"
        )
        # Reads with start (position [0]) >= 100: both reads
        assert count == 2

    def test_fiveprime_plus_strand(self, mock_read_intervals_dict):
        """Test 5' mapping on plus strand (uses read start)"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 150, 250, "fiveprime"
        )
        # Reads with start (position [0]) >= 150: (150,179) and (200,229)
        assert count == 2

    def test_fiveprime_minus_strand(self, mock_read_intervals_dict):
        """Test 5' mapping on minus strand (uses read end)"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "-", 100, 130, "fiveprime"
        )
        # Reads with end (position [1]) <= 130: (100,129)
        assert count == 1

    def test_centered_mapping(self, mock_read_intervals_dict):
        """Test centered mapping uses read center"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 110, 170, "centered"
        )

        assert count == 2

    def test_no_overlapping_reads(self, mock_read_intervals_dict):
        """Test region with no overlapping reads"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 500, 600, "global"
        )
        assert count == 0

    def test_exact_boundary(self, mock_read_intervals_dict):
        """Test reads exactly at boundaries"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 100, 129, "global"
        )
        # First read (100, 129) should overlap
        assert count == 1

    def test_invalid_mapping_method(self, mock_read_intervals_dict):
        """Test invalid mapping method raises no error but counts zero"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 100, 250, "invalid_method"
        )
        assert count == 0

    def test_returns_integer(self, mock_read_intervals_dict):
        """Test that result is an integer"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 100, 250, "global"
        )
        assert isinstance(count, int)

    def test_one_nucleotide_overlap(self, mock_read_intervals_dict):
        """Test reads that overlap by only one nucleotide"""
        count = count_reads(
            mock_read_intervals_dict, "chr1", "+", 129, 130, "global"
        )
        # The read (100,129) overlaps at position 129
        assert count == 1



class TestWindowNormalizeDf:
    """Test suite for window_normalize_df function"""

    @pytest.fixture
    def sample_df(self):
        """Create sample DataFrame for normalization"""
        data = {
            'coordinates': [-50, -25, 0, 25, 50],
            '28': [10, 20, 30, 20, 10],
            '29': [5, 10, 15, 10, 5],
            '30': [20, 40, 60, 40, 20]
        }
        return pd.DataFrame(data)

    def test_normalizes_all_columns_except_first(self, sample_df):
        """Test that all columns except coordinates are normalized"""
        window_size = 100
        result = window_normalize_df(sample_df.copy(), window_size)

        # Coordinates column should be unchanged
        pd.testing.assert_series_equal(result['coordinates'], sample_df['coordinates'])

        # Other columns should be normalized
        assert not result['28'].equals(sample_df['28'])
        assert not result['29'].equals(sample_df['29'])
        assert not result['30'].equals(sample_df['30'])

    def test_normalization_formula(self, sample_df):
        """Test that normalization formula is correct"""
        window_size = 100
        result = window_normalize_df(sample_df.copy(), window_size)

        # For column '28': sum = 90, normalization factor = 90/100 = 0.9
        # Each value should be divided by 0.9
        expected_28 = sample_df['28'] / (sample_df['28'].sum() / window_size)
        pd.testing.assert_series_equal(result['28'], expected_28)

    def test_window_size_affects_normalization(self, sample_df):
        """Test that different window sizes produce different results"""
        result_100 = window_normalize_df(sample_df.copy(), 100)
        result_50 = window_normalize_df(sample_df.copy(), 50)

        # Different window sizes should give different results
        assert not result_100['28'].equals(result_50['28'])

    def test_empty_dataframe(self):
        """Test with empty DataFrame"""
        df = pd.DataFrame({'coordinates': []})
        result = window_normalize_df(df, 100)
        assert len(result) == 0

    def test_single_row(self):
        """Test with single row DataFrame"""
        df = pd.DataFrame({
            'coordinates': [0],
            '30': [100]
        })
        result = window_normalize_df(df, 50)
        # 100 / (100/50) = 100 / 2 = 50
        assert result['30'].iloc[0] == 50.0


class TestJsonDictToDataframe:
    """Test suite for json_dict_to_dataframe function"""

    @pytest.fixture
    def coverage_dict(self):
        """Create sample coverage dictionary"""
        return {
            'chr1': {
                'sample1': {
                    30: {
                        0: {'gene1': 5, 'gene2': 3},
                        25: {'gene1': 2}
                    },
                    32: {
                        0: {'gene1': 10}
                    }
                },
                'sample2': {
                    30: {
                        0: {'gene3': 7}
                    }
                }
            },
            'chr2': {
                'sample1': {
                    30: {
                        -10: {'gene4': 1}
                    }
                }
            }
        }

    def test_creates_correct_columns(self, coverage_dict):
        """Test that DataFrame has correct columns"""
        result = json_dict_to_dataframe(coverage_dict)

        expected_columns = ['chrom', 'sample', 'read_length', 'position', 'gene_id', 'count']
        assert list(result.columns) == expected_columns

    def test_flattens_nested_structure(self, coverage_dict):
        """Test that nested structure is properly flattened"""
        result = json_dict_to_dataframe(coverage_dict)

        # Should have one row per (chrom, sample, read_length, position, gene_id)
        # chr1/sample1/30/0/gene1, chr1/sample1/30/0/gene2, chr1/sample1/30/25/gene1,
        # chr1/sample1/32/0/gene1, chr1/sample2/30/0/gene3, chr2/sample1/30/-10/gene4
        assert len(result) == 6

    def test_converts_to_integers(self, coverage_dict):
        """Test that read_length and position are converted to integers"""
        result = json_dict_to_dataframe(coverage_dict)

        assert result['read_length'].dtype == np.int64
        assert result['position'].dtype == np.int64

    def test_preserves_values(self, coverage_dict):
        """Test that values are correctly preserved"""
        result = json_dict_to_dataframe(coverage_dict)

        # Find row for chr1/sample1/30/0/gene1
        row = result[
            (result['chrom'] == 'chr1') &
            (result['sample'] == 'sample1') &
            (result['read_length'] == 30) &
            (result['position'] == 0) &
            (result['gene_id'] == 'gene1')
        ]

        assert len(row) == 1
        assert row['count'].iloc[0] == 5

    def test_empty_dict(self):
        """Test with empty dictionary"""
        result = json_dict_to_dataframe({})
        assert len(result) == 0
        assert list(result.columns) == ['chrom', 'sample', 'read_length', 'position', 'gene_id', 'count']

class TestAggregateCoverageDataframe:
    """Test suite for aggregate_coverage_dataframe function"""

    @pytest.fixture
    def coverage_df(self):
        """Create comprehensive coverage DataFrame with multiple chromosomes, samples, read lengths, and positions"""
        data = {
            'chrom': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
            'sample': ['sample1', 'sample1', 'sample1', 'sample1', 'sample2', 'sample2', 'sample1', 'sample1', 'sample2'],
            'read_length': [30, 30, 30, 32, 30, 30, 30, 30, 30],
            'position': [0, 0, 5, 0, 0, 5, 0, 5, 0],
            'gene_id': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6', 'gene7', 'gene8', 'gene9'],
            'count': [5, 3, 7, 10, 15, 20, 2, 4, 6]
        }
        return pd.DataFrame(data)

    def test_aggregates_across_genes(self, coverage_df):
        """Test that counts are summed across genes at same position"""
        result = aggregate_coverage_dataframe(coverage_df)

        # chr1/sample1/read_length=30/position=0 should sum gene1(5) + gene2(3) = 8
        row = result[
            (result['chrom'] == 'chr1') &
            (result['sample'] == 'sample1') &
            (result['read_length'] == 30) &
            (result['position'] == 0)
        ]

        assert len(row) == 1
        assert row['count'].iloc[0] == 8

    def test_keeps_separate_positions(self, coverage_df):
        """Test that different positions are kept separate (same chrom/sample/read_length)"""
        result = aggregate_coverage_dataframe(coverage_df)

        # chr1/sample1/read_length=30 should have 2 separate rows for position 0 and 5
        chr1_sample1_30 = result[
            (result['chrom'] == 'chr1') &
            (result['sample'] == 'sample1') &
            (result['read_length'] == 30)
        ]

        assert len(chr1_sample1_30) == 2  # Two positions

        # Position 0: gene1(5) + gene2(3) = 8
        pos_0 = chr1_sample1_30[chr1_sample1_30['position'] == 0]
        assert len(pos_0) == 1
        assert pos_0['count'].iloc[0] == 8

        # Position 5: gene3(7)
        pos_5 = chr1_sample1_30[chr1_sample1_30['position'] == 5]
        assert len(pos_5) == 1
        assert pos_5['count'].iloc[0] == 7

    def test_keeps_separate_chromosomes(self, coverage_df):
        """Test that different chromosomes are kept separate"""
        result = aggregate_coverage_dataframe(coverage_df)

        # Should have rows for both chr1 and chr2
        assert 'chr1' in result['chrom'].values
        assert 'chr2' in result['chrom'].values

        # chr1 should have more rows than chr2
        chr1_rows = result[result['chrom'] == 'chr1']
        chr2_rows = result[result['chrom'] == 'chr2']

        assert len(chr1_rows) == 5
        assert len(chr2_rows) == 3

    def test_keeps_separate_samples(self, coverage_df):
        """Test that different samples are kept separate"""
        result = aggregate_coverage_dataframe(coverage_df)

        # chr1/position=0/read_length=30 should have separate rows for sample1 and sample2
        chr1_pos0_30 = result[
            (result['chrom'] == 'chr1') &
            (result['position'] == 0) &
            (result['read_length'] == 30)
        ]

        # Should have 2 rows: one for sample1, one for sample2
        assert len(chr1_pos0_30) == 2

        # sample1: gene1(5) + gene2(3) = 8
        sample1_row = chr1_pos0_30[chr1_pos0_30['sample'] == 'sample1']
        assert len(sample1_row) == 1
        assert sample1_row['count'].iloc[0] == 8

        # sample2: gene5(15)
        sample2_row = chr1_pos0_30[chr1_pos0_30['sample'] == 'sample2']
        assert len(sample2_row) == 1
        assert sample2_row['count'].iloc[0] == 15

    def test_keeps_separate_read_lengths(self, coverage_df):
        """Test that different read lengths are kept separate"""
        result = aggregate_coverage_dataframe(coverage_df)

        # chr1/sample1/position=0 should have separate rows for read_length 30 and 32
        chr1_sample1_pos0 = result[
            (result['chrom'] == 'chr1') &
            (result['sample'] == 'sample1') &
            (result['position'] == 0)
        ]

        # Should have 2 rows: one for read_length 30, one for read_length 32
        assert len(chr1_sample1_pos0) == 2

        # read_length 30: gene1(5) + gene2(3) = 8
        rl_30 = chr1_sample1_pos0[chr1_sample1_pos0['read_length'] == 30]
        assert len(rl_30) == 1
        assert rl_30['count'].iloc[0] == 8

        # read_length 32: gene4(10)
        rl_32 = chr1_sample1_pos0[chr1_sample1_pos0['read_length'] == 32]
        assert len(rl_32) == 1
        assert rl_32['count'].iloc[0] == 10

    def test_filters_by_read_length(self, coverage_df):
        """Test filtering by read length list"""
        result = aggregate_coverage_dataframe(coverage_df, read_length_list=[30])

        # Should only include read_length 30
        assert all(result['read_length'] == 30)
        assert 32 not in result['read_length'].values

    def test_no_gene_id_column_in_output(self, coverage_df):
        """Test that gene_id column is removed after aggregation"""
        result = aggregate_coverage_dataframe(coverage_df)
        assert 'gene_id' not in result.columns

    def test_output_columns(self, coverage_df):
        """Test that output has correct columns"""
        result = aggregate_coverage_dataframe(coverage_df)
        expected_columns = ['chrom', 'sample', 'read_length', 'position', 'count']
        assert list(result.columns) == expected_columns


    def test_aggregated_config_fills_missing_positions(self, coverage_df):
        """Test that missing positions are filled with zeros based on config"""
        config = {
            'positionsOutsideORF': 10,
            'positionsInsideORF': 10
        }
        result = aggregate_coverage_dataframe(coverage_df, config=config)

        # For chr1/sample1/read_length=30, positions should range from -10 to 10
        chr1_sample1_30 = result[
            (result['chrom'] == 'chr1') &
            (result['sample'] == 'sample1') &
            (result['read_length'] == 30)
        ]

        expected_positions = set(range(-10, 11))
        actual_positions = set(chr1_sample1_30['position'].values)

        assert expected_positions == actual_positions

        # Positions that were not in original data should have count 0
        for pos in expected_positions:
            row = chr1_sample1_30[chr1_sample1_30['position'] == pos]
            if pos not in coverage_df[
                (coverage_df['chrom'] == 'chr1') &
                (coverage_df['sample'] == 'sample1') &
                (coverage_df['read_length'] == 30)
            ]['position'].values:
                assert row['count'].iloc[0] == 0

    def test_empty_dataframe(self):
        """Test with empty DataFrame"""
        df = pd.DataFrame({
            'chrom': [],
            'sample': [],
            'read_length': [],
            'position': [],
            'gene_id': [],
            'count': []
        })

        result = aggregate_coverage_dataframe(df)

        # Should return empty DataFrame with correct columns
        assert len(result) == 0
        assert list(result.columns) == ['chrom', 'sample', 'read_length', 'position', 'count']

    def test_empty_after_filtering(self, coverage_df):
        """Test when filtering removes all rows"""
        result = aggregate_coverage_dataframe(coverage_df, read_length_list=[99])

        # No reads with length 99, should be empty
        assert len(result) == 0
        assert list(result.columns) == ['chrom', 'sample', 'read_length', 'position', 'count']


class TestApplyNormalization:
    """Test suite for apply_normalization function"""

    @pytest.fixture
    def sample_df(self):
        """Create sample DataFrame for normalization"""
        data = {
            'chrom': ['chr1', 'chr1', 'chr1', 'chr1'],
            'sample': ['sample1', 'sample1', 'sample1', 'sample1'],
            'read_length': [30, 30, 30, 30],
            'position': [0, 25, 50, 75],
            'count': [100, 200, 300, 400]
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def total_counts_dict(self):
        """Create sample total counts dictionary"""
        return {
            'sample1': {'chr1': 1000000},
            'sample2': {'chr2': 500000}
        }

    def test_raw_normalization_unchanged(self, sample_df):
        """Test that raw normalization returns unchanged data"""
        result = apply_normalization(sample_df, "raw")
        pd.testing.assert_frame_equal(result, sample_df)

    def test_cpm_normalization_formula(self, sample_df, total_counts_dict):
        """Test CPM normalization formula"""
        result = apply_normalization(sample_df, "cpm", total_counts_dict)

        # CPM = count / total_counts * 1e6
        # For chr1/sample1: total = 1000000
        # 100 / 1000000 * 1e6 = 100
        expected_first = 100 / 1000000 * 1e6
        assert result['count'].iloc[0] == expected_first

    def test_cpm_requires_total_counts(self, sample_df):
        """Test that CPM raises error without total_counts_dict"""
        with pytest.raises(ValueError, match="CPM normalization requires"):
            apply_normalization(sample_df, "cpm")

    def test_window_normalization(self, sample_df):
        """Test window normalization"""
        result = apply_normalization(sample_df, "window")

        # Total = 100+200+300+400 = 1000
        # Window size = 4 positions
        # Normalization factor = 1000 / 4 = 250
        # First value: 100 / 250 = 0.4
        expected_first = 100 / (1000 / 4)
        assert result['count'].iloc[0] == expected_first

    def test_returns_copy(self, sample_df):
        """Test that normalization returns a copy, not modifying original"""
        original_values = sample_df['count'].copy()
        result = apply_normalization(sample_df, "cpm", {'sample1': 1000000})

        # Original should be unchanged
        pd.testing.assert_series_equal(sample_df['count'], original_values)

    def test_cpm_with_simple_total_counts(self, sample_df):
        """Test CPM with total counts not per chromosome"""
        simple_counts = {'sample1': 1000000}
        result = apply_normalization(sample_df, "cpm", simple_counts)

        assert len(result) == len(sample_df)

    def test_window_normalization_groups_correctly(self):
        """Test that window normalization groups by chrom/sample/read_length"""
        data = {
            'chrom': ['chr1', 'chr1', 'chr2', 'chr2'],
            'sample': ['sample1', 'sample1', 'sample1', 'sample1'],
            'read_length': [30, 30, 30, 30],
            'position': [0, 25, 0, 25],
            'count': [100, 100, 50, 50]
        }
        df = pd.DataFrame(data)

        result = apply_normalization(df, "window")

        # Each chromosome should be normalized separately
        # chr1: total=200, window=2, factor=100, values = 100/100 = 1.0
        chr1_values = result[result['chrom'] == 'chr1']['count'].values
        assert all(chr1_values == 1.0)

    def test_invalid_normalization_method(self, sample_df):
        """Test that invalid normalization method raises no error but returns unchanged"""
        result = apply_normalization(sample_df, "invalid_method")
        pd.testing.assert_frame_equal(result, sample_df)

    def test_empty_dataframe(self):
        """Test with empty DataFrame"""
        df = pd.DataFrame({
            'chrom': [],
            'sample': [],
            'read_length': [],
            'position': [],
            'count': []
        })

        result = apply_normalization(df, "cpm", {'sample1': 1000000})

        # Should return empty DataFrame with correct columns
        assert len(result) == 0
        assert list(result.columns) == ['chrom', 'sample', 'read_length', 'position', 'count']

if __name__ == "__main__":
    pytest.main([__file__, "-v"])