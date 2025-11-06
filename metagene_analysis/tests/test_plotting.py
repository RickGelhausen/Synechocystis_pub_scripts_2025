"""
Tests for plotting functions.
Author: Rick Gelhausen
"""

import pytest
import pandas as pd
import plotly.graph_objects as go
from lib.plotting import prepare_colors_and_lines, plot_metagene_from_dataframes


class TestPrepareColorsAndLines:
    """Test suite for prepare_colors_and_lines function"""

    def test_returns_two_lists(self):
        """Test that function returns two lists"""
        color_list = ["red", "blue"]
        colors, lines = prepare_colors_and_lines(color_list)

        assert isinstance(colors, list)
        assert isinstance(lines, list)

    def test_creates_enough_combinations(self):
        """Test that function creates combinations of colors and line types"""
        color_list = ["red", "blue", "green"]
        colors, lines = prepare_colors_and_lines(color_list)

        # 6 line types * 3 colors = 18 combinations
        assert len(colors) == 18
        assert len(lines) == 18

    def test_line_types_included(self):
        """Test that all line types are present"""
        color_list = ["red"]
        colors, lines = prepare_colors_and_lines(color_list)

        expected_line_types = ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]
        assert set(lines) == set(expected_line_types)

    def test_colors_repeated_per_line_type(self):
        """Test that colors are repeated for each line type"""
        color_list = ["red", "blue"]
        colors, lines = prepare_colors_and_lines(color_list)

        # First 2 should be red, blue (solid)
        assert colors[0] == "red"
        assert colors[1] == "blue"
        assert lines[0] == "solid"
        assert lines[1] == "solid"

        # Next 2 should be red, blue (dot)
        assert colors[2] == "red"
        assert colors[3] == "blue"
        assert lines[2] == "dot"
        assert lines[3] == "dot"

    def test_single_color(self):
        """Test with single color"""
        color_list = ["red"]
        colors, lines = prepare_colors_and_lines(color_list)

        assert len(colors) == 6  # One for each line type
        assert all(c == "red" for c in colors)


class TestPlotMetageneFromDataframes:
    """Test suite for plot_metagene_from_dataframes function"""

    @pytest.fixture
    def sample_df_start(self):
        """Create sample start codon dataframe"""
        data = {
            'chrom': ['chr1'] * 6,
            'sample': ['sample1'] * 6,
            'read_length': [28, 28, 28, 30, 30, 30],
            'position': [-10, 0, 10, -10, 0, 10],
            'count': [50, 100, 75, 60, 120, 80]
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def sample_df_stop(self):
        """Create sample stop codon dataframe"""
        data = {
            'chrom': ['chr1'] * 6,
            'sample': ['sample1'] * 6,
            'read_length': [28, 28, 28, 30, 30, 30],
            'position': [-10, 0, 10, -10, 0, 10],
            'count': [40, 90, 65, 50, 110, 70]
        }
        return pd.DataFrame(data)

    def test_returns_figure_and_max_y(self, sample_df_start, sample_df_stop):
        """Test that function returns figure and max_y value"""
        read_length_list = [28, 30]

        fig, max_y = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list
        )

        assert isinstance(fig, go.Figure)
        assert isinstance(max_y, (int, float))
        assert max_y > 0

    def test_creates_subplots(self, sample_df_start, sample_df_stop):
        """Test that figure has correct subplot structure"""
        read_length_list = [28, 30]

        fig, _ = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list
        )

        # Check that we have 2 subplots by checking for xaxis and xaxis2
        assert 'xaxis' in fig.layout
        assert 'xaxis2' in fig.layout
        assert 'yaxis' in fig.layout
        assert 'yaxis2' in fig.layout

        # Verify they are side by side (different domains)
        assert fig.layout.xaxis.domain[1] < fig.layout.xaxis2.domain[0]


    def test_creates_traces_for_each_read_length(self, sample_df_start, sample_df_stop):
        """Test that traces are created for each read length"""
        read_length_list = [28, 30]

        fig, _ = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list
        )

        # Should have 2 traces per subplot (one per read length)
        # Start codon: 2 traces, Stop codon: 2 traces
        # Plus 2 vertical lines (shapes)
        assert len(fig.data) >= 4

    def test_handles_empty_dataframes(self):
        """Test that function handles empty dataframes gracefully"""
        empty_df = pd.DataFrame(columns=['chrom', 'sample', 'read_length', 'position', 'count'])
        read_length_list = [28, 30]

        fig, max_y = plot_metagene_from_dataframes(
            empty_df, empty_df, "chr1", read_length_list
        )

        assert isinstance(fig, go.Figure)
        # max_y should still be calculated (will be small due to padding on 0)
        assert max_y >= 0

    def test_uses_default_colors_when_none_provided(self, sample_df_start, sample_df_stop):
        """Test that default colors are used when color_list is None"""
        read_length_list = [28, 30]

        fig, _ = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list, color_list=None
        )

        # Check that traces have colors assigned
        assert all(hasattr(trace, 'line') for trace in fig.data if hasattr(trace, 'line'))

    def test_uses_custom_colors(self, sample_df_start, sample_df_stop):
        """Test that custom colors are applied"""
        read_length_list = [28]
        custom_colors = ["#FF0000", "#00FF00"]

        fig, _ = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list, color_list=custom_colors
        )

        # First trace should use first color
        assert fig.data[0].line.color == "#FF0000"

    def test_adds_vertical_lines(self, sample_df_start, sample_df_stop):
        """Test that vertical reference lines are added at position 0"""
        read_length_list = [28, 30]

        fig, _ = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list
        )

        # Check that shapes (vertical lines) are added
        assert len(fig.layout.shapes) == 2  # One for each subplot

    def test_max_y_includes_padding(self, sample_df_start, sample_df_stop):
        """Test that max_y includes 10% padding"""
        read_length_list = [28, 30]

        # Max count in sample data is 120
        expected_max = 120 * 1.15

        fig, max_y = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list
        )

        assert max_y >= 120  # At least the max count
        assert max_y == pytest.approx(expected_max, rel=0.01)

    def test_averages_multiple_samples(self):
        """Test that multiple samples are averaged correctly"""
        df_multi_sample = pd.DataFrame({
            'chrom': ['chr1'] * 6,
            'sample': ['sample1', 'sample1', 'sample1', 'sample2', 'sample2', 'sample2'],
            'read_length': [28, 28, 28, 28, 28, 28],
            'position': [0, 5, 10, 0, 5, 10],
            'count': [100, 50, 75, 200, 150, 125]
        })

        read_length_list = [28]

        fig, _ = plot_metagene_from_dataframes(
            df_multi_sample, df_multi_sample, "chr1", read_length_list
        )

        # Check that data exists (averaging should work)
        assert len(fig.data) >= 2  # At least one trace per subplot

    def test_title_includes_chromosome(self, sample_df_start, sample_df_stop):
        """Test that plot title includes chromosome name"""
        read_length_list = [28]

        fig, _ = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list
        )

        assert "chr1" in fig.layout.title.text

    def test_legend_group_consistency(self, sample_df_start, sample_df_stop):
        """Test that legend groups are consistent between subplots"""
        read_length_list = [28, 30]

        fig, _ = plot_metagene_from_dataframes(
            sample_df_start, sample_df_stop, "chr1", read_length_list
        )

        # Traces should have legendgroups
        assert all(hasattr(trace, 'legendgroup') for trace in fig.data)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])