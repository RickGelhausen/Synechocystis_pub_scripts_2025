"""
Unit tests for lib.io module

Tests for:
- write_plots_to_file
- create_interactive_html
- parse_bam_file_names
- parse_read_lengths
- parse_genome_lengths
- save_dataframes_to_excel
- excel_writer
"""

import tempfile
import os
import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import Mock, patch, mock_open


from lib.io import (
    write_plots_to_file,
    create_interactive_html,
    parse_bam_file_names,
    parse_read_lengths,
    parse_genome_lengths,
    save_dataframes_to_excel,
    excel_writer,
    INTRO_HTML,
    OUTRO_HTML
)


class TestWritePlotsToFile:
    """Test suite for write_plots_to_file function"""

    @pytest.fixture
    def mock_fig_list(self):
        """Create mock figure list"""
        fig1 = Mock()
        fig2 = Mock()
        fig3 = Mock()
        return [
            ("chr1", "global", fig1),
            ("chr1", "threeprime", fig2),
            ("chr2", "fiveprime", fig3),
        ]

    def test_writes_png_files(self, mock_fig_list):
        """Test PNG output format"""
        output_format = ["png"]

        write_plots_to_file(
            mock_fig_list, output_format, "integrated",
            "test_sample", "/output/meta", 1400, 600
        )

        # Verify write_image was called for each figure
        for chrom, method, fig in mock_fig_list:
            fig.write_image.assert_called_once_with(
                f"/output/meta/{chrom}_{method}.png",
                width=1400,
                height=600
            )

    def test_writes_jpg_files(self, mock_fig_list):
        """Test JPG output format"""
        output_format = ["jpg"]

        write_plots_to_file(
            mock_fig_list, output_format, "integrated",
            "test_sample", "/output/meta"
        )

        for chrom, method, fig in mock_fig_list:
            fig.write_image.assert_called_once_with(
                f"/output/meta/{chrom}_{method}.jpg",
                width=1400,
                height=600
            )

    def test_writes_svg_files(self, mock_fig_list):
        """Test SVG output format"""
        output_format = ["svg"]

        write_plots_to_file(
            mock_fig_list, output_format, "integrated",
            "test_sample", "/output/meta"
        )

        for chrom, method, fig in mock_fig_list:
            fig.write_image.assert_called_once_with(
                f"/output/meta/{chrom}_{method}.svg",
                width=1400,
                height=600
            )

    def test_writes_multiple_formats(self, mock_fig_list):
        """Test multiple output formats"""
        output_format = ["png", "svg"]

        write_plots_to_file(
            mock_fig_list, output_format, "integrated",
            "test_sample", "/output/meta"
        )

        # Each figure should be written twice (png and svg)
        for _, _, fig in mock_fig_list:
            assert fig.write_image.call_count == 2

    @patch('lib.io.create_interactive_html')
    def test_creates_interactive_html(self, mock_create_html, mock_fig_list):
        """Test interactive HTML output"""
        output_format = ["interactive"]

        write_plots_to_file(
            mock_fig_list, output_format, "integrated",
            "test_sample", "/output/meta"
        )

        mock_create_html.assert_called_once_with(
            mock_fig_list,
            "test_sample",
            "/output/meta/interactive_metagene_profiling.html",
            "integrated"
        )

    def test_custom_dimensions(self, mock_fig_list):
        """Test custom figure dimensions"""
        output_format = ["png"]

        write_plots_to_file(
            mock_fig_list, output_format, "integrated",
            "test_sample", "/output/meta",
            fig_width=2000, fig_height=800
        )

        mock_fig_list[0][2].write_image.assert_called_once_with(
            "/output/meta/chr1_global.png",
            width=2000,
            height=800
        )

    def test_no_output_format(self, mock_fig_list):
        """Test with empty output format list"""
        output_format = []

        # Should not raise error, just do nothing
        write_plots_to_file(
            mock_fig_list, output_format, "integrated",
            "test_sample", "/output/meta"
        )

        # No figures should be written
        for _, _, fig in mock_fig_list:
            fig.write_image.assert_not_called()


class TestCreateInteractiveHtml:
    """Test suite for create_interactive_html function"""

    @pytest.fixture
    def mock_fig_list(self):
        """Create mock figure list with different mapping methods"""
        fig1 = Mock()
        fig1.to_html.return_value = "<div>plot1</div>"
        fig2 = Mock()
        fig2.to_html.return_value = "<div>plot2</div>"
        fig3 = Mock()
        fig3.to_html.return_value = "<div>plot3</div>"

        return [
            ("chr1", "global", fig1),
            ("chr1", "threeprime", fig2),
            ("chr2", "fiveprime", fig3),
        ]

    @patch('builtins.open', new_callable=mock_open)
    def test_creates_html_file(self, mock_file, mock_fig_list):
        """Test that HTML file is created"""
        create_interactive_html(
            mock_fig_list, "test_sample", "/output/test.html", "integrated"
        )

        mock_file.assert_called_once_with("/output/test.html", "w", encoding="utf-8")
        mock_file().write.assert_called_once()

    @patch('builtins.open', new_callable=mock_open)
    def test_includes_intro_and_outro(self, mock_file, mock_fig_list):
        """Test that INTRO_HTML and OUTRO_HTML are included"""
        create_interactive_html(
            mock_fig_list, "test_sample", "/output/test.html", "integrated"
        )

        written_content = mock_file().write.call_args[0][0]
        assert INTRO_HTML in written_content
        assert OUTRO_HTML in written_content

    @patch('builtins.open', new_callable=mock_open)
    def test_includes_alignment_name(self, mock_file, mock_fig_list):
        """Test that alignment file name is included as header"""
        create_interactive_html(
            mock_fig_list, "my_sample_name", "/output/test.html", "integrated"
        )

        written_content = mock_file().write.call_args[0][0]
        assert "<h1>my_sample_name</h1>" in written_content

    @patch('builtins.open', new_callable=mock_open)
    def test_mapping_method_headers(self, mock_file, mock_fig_list):
        """Test that mapping method headers are created"""
        create_interactive_html(
            mock_fig_list, "test", "/output/test.html", "integrated"
        )

        written_content = mock_file().write.call_args[0][0]
        assert "<h2>Global Mapping</h2>" in written_content
        assert "<h2>3' Mapping</h2>" in written_content
        assert "<h2>5' Mapping</h2>" in written_content

    @patch('builtins.open', new_callable=mock_open)
    def test_centered_mapping_label(self, mock_file):
        """Test centered mapping method label"""
        fig = Mock()
        fig.to_html.return_value = "<div>plot</div>"
        fig_list = [("chr1", "centered", fig)]

        create_interactive_html(fig_list, "test", "/output/test.html", "integrated")

        written_content = mock_file().write.call_args[0][0]
        assert "<h2>Centered Mapping</h2>" in written_content

    @patch('builtins.open', new_callable=mock_open)
    def test_unknown_mapping_label(self, mock_file):
        """Test unknown mapping method gets labeled as Unknown"""
        fig = Mock()
        fig.to_html.return_value = "<div>plot</div>"
        fig_list = [("chr1", "weird_method", fig)]

        create_interactive_html(fig_list, "test", "/output/test.html", "integrated")

        written_content = mock_file().write.call_args[0][0]
        assert "<h2>Unknown</h2>" in written_content

    @patch('builtins.open', new_callable=mock_open)
    def test_integrated_plotly_js(self, mock_file, mock_fig_list):
        """Test integrated plotly.js mode"""
        create_interactive_html(
            mock_fig_list, "test", "/output/test.html", "integrated"
        )

        # First figure should include plotly.js
        mock_fig_list[0][2].to_html.assert_called_once()
        call_kwargs = mock_fig_list[0][2].to_html.call_args[1]
        assert call_kwargs['full_html'] is False
        assert 'include_plotlyjs' not in call_kwargs  # Default behavior

        # Other figures should exclude plotly.js
        mock_fig_list[1][2].to_html.assert_called_once()
        call_kwargs = mock_fig_list[1][2].to_html.call_args[1]
        assert call_kwargs['include_plotlyjs'] is False

    @patch('builtins.open', new_callable=mock_open)
    def test_online_plotly_js(self, mock_file, mock_fig_list):
        """Test online CDN plotly.js mode"""
        create_interactive_html(
            mock_fig_list, "test", "/output/test.html", "online"
        )

        # First figure should use CDN
        call_kwargs = mock_fig_list[0][2].to_html.call_args[1]
        assert call_kwargs['include_plotlyjs'] == "cdn"

    @patch('builtins.open', new_callable=mock_open)
    def test_local_plotly_js(self, mock_file, mock_fig_list):
        """Test local directory plotly.js mode"""
        create_interactive_html(
            mock_fig_list, "test", "/output/test.html", "local"
        )

        # First figure should use directory
        call_kwargs = mock_fig_list[0][2].to_html.call_args[1]
        assert call_kwargs['include_plotlyjs'] == "directory"

    @patch('builtins.open', new_callable=mock_open)
    def test_svg_download_config(self, mock_file, mock_fig_list):
        """Test that SVG download is configured"""
        create_interactive_html(
            mock_fig_list, "test", "/output/test.html", "integrated"
        )

        # Check that config includes SVG format
        for _, _, fig in mock_fig_list:
            call_kwargs = fig.to_html.call_args[1]
            assert call_kwargs['config']['toImageButtonOptions']['format'] == 'svg'


class TestParseBamFileNames:
    """Test suite for parse_bam_file_names function"""

    @pytest.fixture
    def mock_bam_folder(self):
        """Create mock folder with BAM files"""
        folder = Mock(spec=Path)

        # Create mock BAM files
        bam1 = Mock(spec=Path)
        bam1.name = "sample1.bam"

        bam2 = Mock(spec=Path)
        bam2.name = "sample2.bam"

        bam3 = Mock(spec=Path)
        bam3.name = "sample3.sorted.bam"

        folder.glob.return_value = [bam1, bam2, bam3]

        return folder

    def test_returns_dict_with_bam_files(self, mock_bam_folder):
        """Test that function returns dictionary of BAM files"""
        result = parse_bam_file_names(mock_bam_folder)

        assert isinstance(result, dict)
        assert len(result) == 3

    def test_extracts_prefix_correctly(self, mock_bam_folder):
        """Test that prefix (before first dot) is used as key"""
        result = parse_bam_file_names(mock_bam_folder)

        assert "sample1" in result
        assert "sample2" in result
        assert "sample3" in result  # Should use "sample3", not "sample3.sorted"

    def test_stores_path_objects(self, mock_bam_folder):
        """Test that Path objects are stored as values"""
        result = parse_bam_file_names(mock_bam_folder)

        # Values should be Path objects (our mocks are spec'd as Path)
        for value in result.values():
            # Check it has Path-like attributes
            assert hasattr(value, 'name')
            assert hasattr(value, 'parent')

    def test_glob_called_correctly(self, mock_bam_folder):
        """Test that glob is called with correct pattern"""
        parse_bam_file_names(mock_bam_folder)

        mock_bam_folder.glob.assert_called_once_with("*.bam")

    def test_empty_folder(self):
        """Test with folder containing no BAM files"""
        folder = Mock(spec=Path)
        folder.glob.return_value = []

        result = parse_bam_file_names(folder)

        assert not result


class TestParseReadLengths:
    """Test suite for parse_read_lengths function"""

    def test_single_value(self):
        """Test parsing single read length"""
        result = parse_read_lengths("30")
        assert result == [30]

    def test_multiple_values(self):
        """Test parsing comma-separated values"""
        result = parse_read_lengths("28,29,30,31")
        assert result == [28, 29, 30, 31]

    def test_range_ascending(self):
        """Test parsing ascending range"""
        result = parse_read_lengths("28-32")
        assert result == [28, 29, 30, 31, 32]

    def test_range_descending(self):
        """Test parsing descending range (should be reversed)"""
        result = parse_read_lengths("32-28")
        assert result == [28, 29, 30, 31, 32]

    def test_mixed_values_and_ranges(self):
        """Test parsing mix of values and ranges"""
        result = parse_read_lengths("25,28-30,35")
        assert result == [25, 28, 29, 30, 35]

    def test_removes_duplicates(self):
        """Test that duplicate values are removed"""
        result = parse_read_lengths("28,29,28-30")
        assert result == [28, 29, 30]

    def test_sorted_output(self):
        """Test that output is sorted"""
        result = parse_read_lengths("30,25,28")
        assert result == [25, 28, 30]

    def test_single_element_range(self):
        """Test range with same start and end"""
        result = parse_read_lengths("30-30")
        assert result == [30]

    def test_complex_input(self):
        """Test complex input with multiple ranges and values"""
        result = parse_read_lengths("20-22,25,28-30,32,35-37")
        assert result == [20, 21, 22, 25, 28, 29, 30, 32, 35, 36, 37]

    def test_returns_list_not_set(self):
        """Test that return type is list"""
        result = parse_read_lengths("28,29,30")
        assert isinstance(result, list)

    def test_empty_string(self):
        """Test with empty input string"""

        with pytest.raises(ValueError, match="Read lengths must be numeric values."):
            parse_read_lengths("a")


class TestParseGenomeLengths:
    """Test suite for parse_genome_lengths function"""

    @pytest.fixture
    def genome_file_content(self):
        """Sample FASTA file content"""
        return """>chr1 chromosome 1
ATCGATCGATCG
GCTAGCTAGCTA
>chr2 chromosome 2
ATCG
GCTA
>chr3
ATCGATCGATCGATCG
"""

    def test_parses_chromosome_names(self, genome_file_content):
        """Test that chromosome names are extracted correctly"""
        with patch('lib.io.open', mock_open(read_data=genome_file_content)):
            result = parse_genome_lengths("genome.fa")
        assert "chr1" in result
        assert "chr2" in result
        assert "chr3" in result

    def test_calculates_lengths_correctly(self, genome_file_content):
        """Test that sequence lengths are calculated correctly"""
        with patch('lib.io.open', mock_open(read_data=genome_file_content)):
            result = parse_genome_lengths("genome.fa")
        # chr1: 12 + 12 = 24
        assert result["chr1"] == 24
        # chr2: 4 + 4 = 8
        assert result["chr2"] == 8
        # chr3: 16
        assert result["chr3"] == 16

    def test_ignores_header_annotations(self, genome_file_content):
        """Test that text after chromosome name in header is ignored"""
        with patch('lib.io.open', mock_open(read_data=genome_file_content)):
            result = parse_genome_lengths("genome.fa")
        # Should be "chr1", not "chr1 chromosome 1"
        assert "chr1" in result
        assert "chr1 chromosome 1" not in result

    def test_strips_whitespace(self):
        """Test that whitespace in sequences is stripped"""
        content = """>chr1
ATCG
  GCTA
"""
        with patch('lib.io.open', mock_open(read_data=content)):
            result = parse_genome_lengths("genome.fa")
        # Should not count the spaces
        assert result["chr1"] == 8

    def test_empty_file(self):
        """Test with empty file"""
        content = ""
        with patch('lib.io.open', mock_open(read_data=content)):
            result = parse_genome_lengths("genome.fa")
        assert not result

    def test_file_encoding(self):
        """Test that file is opened with UTF-8 encoding"""
        content = ">chr1\nATCG\n"
        m = mock_open(read_data=content)
        with patch('lib.io.open', m):
            result = parse_genome_lengths("genome.fa")
        m.assert_called_once_with("genome.fa", "r", encoding="utf-8")


class TestSaveDataframesToExcel:
    """Test suite for save_dataframes_to_excel function"""

    @pytest.fixture
    def sample_df(self):
        """Create sample aggregated DataFrame for a single sample"""
        data = {
            'chrom': ['chr1', 'chr1', 'chr1', 'chr2', 'chr2'],
            'sample': ['sample1', 'sample1', 'sample1', 'sample1', 'sample1'],  # All same sample
            'read_length': [30, 30, 32, 30, 32],  # Different read lengths
            'position': [0, 25, 0, 0, 0],
            'count': [100, 200, 50, 10, 20]
        }
        return pd.DataFrame(data)

    @patch('lib.io.excel_writer')
    def test_calls_excel_writer(self, mock_excel_writer, sample_df):
        """Test that excel_writer is called"""
        save_dataframes_to_excel(sample_df, "output.xlsx")

        mock_excel_writer.assert_called_once()
        assert mock_excel_writer.call_args[0][0] == "output.xlsx"

    @patch('lib.io.excel_writer')
    def test_creates_sheet_per_chromosome(self, mock_excel_writer, sample_df):
        """Test that separate sheet is created for each chromosome"""
        save_dataframes_to_excel(sample_df, "output.xlsx")

        df_dict = mock_excel_writer.call_args[0][1]
        assert 'chr1' in df_dict
        assert 'chr2' in df_dict

    @patch('lib.io.excel_writer')
    def test_pivots_data_correctly(self, mock_excel_writer, sample_df):
        """Test that data is pivoted with position as rows"""
        save_dataframes_to_excel(sample_df, "output.xlsx")

        df_dict = mock_excel_writer.call_args[0][1]
        chr1_df = df_dict['chr1']

        # Should have position column
        assert 'position' in chr1_df.columns
        # Position should be sorted
        assert chr1_df['position'].tolist() == [0, 25]

    @patch('lib.io.excel_writer')
    def test_creates_sum_column(self, mock_excel_writer, sample_df):
        """Test that sum column is added"""
        save_dataframes_to_excel(sample_df, "output.xlsx")

        df_dict = mock_excel_writer.call_args[0][1]
        chr1_df = df_dict['chr1']

        assert 'sum' in chr1_df.columns

    @patch('lib.io.excel_writer')
    def test_sum_calculation(self, mock_excel_writer, sample_df):
        """Test that sum is calculated correctly for single sample"""
        save_dataframes_to_excel(sample_df, "output.xlsx")
        df_dict = mock_excel_writer.call_args[0][1]
        chr1_df = df_dict['chr1']

        # At position 0: read_length columns sum to expected total
        pos_0_row = chr1_df[chr1_df['position'] == 0]
        assert pos_0_row['sum'].iloc[0] == 150

    @patch('lib.io.excel_writer')
    def test_rejects_multiple_samples(self, mock_excel_writer):
        """Test that function rejects multi-sample dataframes"""
        multi_sample_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'position': [0, 0],
            'sample': ['sample1', 'sample2'],
            'read_length': [100, 100],
            'count': [50, 50]
        })

        with pytest.raises(ValueError, match="expects single sample"):
            save_dataframes_to_excel(multi_sample_df, "output.xlsx")

    @patch('lib.io.excel_writer')
    def test_column_names_formatted(self, mock_excel_writer, sample_df):
        """Test that column names are formatted as read lengths"""
        save_dataframes_to_excel(sample_df, "output.xlsx")

        df_dict = mock_excel_writer.call_args[0][1]
        chr1_df = df_dict['chr1']

        # Should have columns like "30" (as strings)
        columns = [col for col in chr1_df.columns if col not in ['position', 'sum']]
        assert all(isinstance(col, str) for col in columns)
        assert all(col.isdigit() for col in columns)


class TestExcelWriter:
    """Test suite for excel_writer function"""

    @pytest.fixture
    def sample_dataframes(self):
        """Create sample DataFrames for Excel output"""
        df1 = pd.DataFrame({
            'position': [0, 25, 50],
            '28': [10, 20, 30],
            '30': [15, 25, 35],
            'sum': [25, 45, 65]
        })
        df2 = pd.DataFrame({
            'position': [0, 10],
            '30': [100, 200],
            'sum': [100, 200]
        })
        return {'chr1': df1, 'chr2': df2}

    def test_creates_excel_file(self, sample_dataframes):
        """Test that Excel file is created with correct structure"""
        with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as tmp:
            tmp_path = tmp.name

        try:
            # Call the function
            excel_writer(tmp_path, sample_dataframes)

            # Verify file was created
            assert os.path.exists(tmp_path)

            # Read back and verify content
            excel_data = pd.read_excel(tmp_path, sheet_name=None)
            assert 'chr1' in excel_data
            assert 'chr2' in excel_data

            # Verify chr1 data
            chr1_df = excel_data['chr1']
            assert len(chr1_df) == 3
            assert list(chr1_df['position']) == [0, 25, 50]

            # Verify chr2 data
            chr2_df = excel_data['chr2']
            assert len(chr2_df) == 2

        finally:
            # Cleanup
            if os.path.exists(tmp_path):
                os.remove(tmp_path)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])