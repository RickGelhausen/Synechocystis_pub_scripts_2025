"""
Integration tests for metagene workflow.
"""

import copy
from pathlib import Path
import tempfile
import pandas as pd
import pytest
import pysam
from metagene_profiling import create_metagene_plots, extract_read_counts

class TestCreateMetagenePlotsIntegration:
    """Integration tests for create_metagene_plots workflow"""

    @pytest.fixture
    def sample_results(self):
        """Create sample results dictionary with correct nested structure"""
        return {
            'threeprime': {
                'start': {
                    'chr1': {
                        'sample1': {
                            28: {
                                0: {'gene1': 100, 'gene2': 50},
                                5: {'gene1': 50, 'gene2': 25},
                                10: {'gene1': 75, 'gene2': 30}
                            },
                            30: {
                                0: {'gene1': 120, 'gene2': 60},
                                5: {'gene1': 60, 'gene2': 30},
                                10: {'gene1': 80, 'gene2': 40}
                            }
                        }
                    },
                    'chr2': {
                        'sample1': {
                            28: {
                                0: {'gene3': 40},
                                5: {'gene3': 20}
                            },
                            30: {
                                0: {'gene3': 50},
                                5: {'gene3': 25}
                            }
                        }
                    }
                },
                'stop': {
                    'chr1': {
                        'sample1': {
                            28: {
                                0: {'gene1': 90, 'gene2': 45},
                                5: {'gene1': 45, 'gene2': 20},
                                10: {'gene1': 65, 'gene2': 25}
                            },
                            30: {
                                0: {'gene1': 110, 'gene2': 55},
                                5: {'gene1': 55, 'gene2': 25},
                                10: {'gene1': 70, 'gene2': 35}
                            }
                        }
                    },
                    'chr2': {
                        'sample1': {
                            28: {
                                0: {'gene3': 35},
                                5: {'gene3': 18}
                            },
                            30: {
                                0: {'gene3': 45},
                                5: {'gene3': 22}
                            }
                        }
                    }
                },
                'total_counts': {
                    'sample1': {28: 10000, 30: 12000}
                }
            }
        }

    @pytest.fixture
    def sample_config(self):
        """Create sample configuration"""
        return {
            'positionsOutsideORF': 100,
            'positionsInsideORF': 250,
            'filteringMethods': [],
            'neighboringGenesDistance': 50,
            'rpkmThreshold': 10.0,
            'lengthCutoff': 50,
            'mappingMethods': ['threeprime'],
            'readLengths': '28-30',
            'colorList': [],
            'normalizationMethods': ['raw'],
            'outputFormats': ['interactive'],
            'includePlotlyJS': "integrated"
        }

    def test_creates_output_directory(self, sample_results, sample_config):
        """Test that output directories are created"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Check that directories were created
            assert output_dir.exists()
            assert (output_dir / "raw").exists()

    def test_creates_excel_files(self, sample_results, sample_config):
        """Test that Excel files are created for each mapping mode"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Check Excel files exist
            assert (output_dir / "raw" / "threeprime_readcounts_start.xlsx").exists()
            assert (output_dir / "raw" / "threeprime_readcounts_stop.xlsx").exists()

    def test_excel_has_correct_sheet_structure(self, sample_results, sample_config):
        """Test that Excel files have one sheet per chromosome"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Read Excel file
            excel_path = output_dir / "raw" / "threeprime_readcounts_start.xlsx"
            df_dict = pd.read_excel(excel_path, sheet_name=None)

            # Should have sheets for chr1 and chr2
            assert 'chr1' in df_dict
            assert 'chr2' in df_dict

    def test_excel_contains_expected_columns(self, sample_results, sample_config):
        """Test that Excel sheets have expected column structure"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Read Excel file
            excel_path = output_dir / "raw" / "threeprime_readcounts_start.xlsx"
            df_dict = pd.read_excel(excel_path, sheet_name=None)
            chr1_df = df_dict['chr1']

            # Should have position column
            assert 'position' in chr1_df.columns

            # Should have sum column
            assert 'sum' in chr1_df.columns

            # Should have read length columns (28, 30)
            assert '28' in chr1_df.columns or 28 in chr1_df.columns
            assert '30' in chr1_df.columns or 30 in chr1_df.columns

    def test_excel_aggregates_across_genes(self, sample_results, sample_config):
        """Test that Excel data aggregates counts across genes"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Read Excel file
            excel_path = output_dir / "raw" / "threeprime_readcounts_start.xlsx"
            df_dict = pd.read_excel(excel_path, sheet_name=None)
            chr1_df = df_dict['chr1']

            # At position 0, read_length 28: gene1(100) + gene2(50) = 150
            # Get row for position 0
            pos_0 = chr1_df[chr1_df['position'] == 0]
            assert not pos_0.empty, "Position 0 should exist in dataframe"

            # Check that read length 28 column has aggregated value
            # (exact value depends on how columns are named)
            if '28' in chr1_df.columns:
                assert pos_0['28'].iloc[0] == 150
            elif 28 in chr1_df.columns:
                assert pos_0[28].iloc[0] == 150

    def test_excel_sum_column_calculation(self, sample_results, sample_config):
        """Test that sum column correctly sums across read lengths"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Read Excel file
            excel_path = output_dir / "raw" / "threeprime_readcounts_start.xlsx"
            df_dict = pd.read_excel(excel_path, sheet_name=None)
            chr1_df = df_dict['chr1']

            # At position 0:
            # read_length 28: 150 (gene1:100 + gene2:50)
            # read_length 30: 180 (gene1:120 + gene2:60)
            # sum should be: 330
            pos_0 = chr1_df[chr1_df['position'] == 0]
            assert pos_0['sum'].iloc[0] == 330

    def test_creates_html_plots(self, sample_results, sample_config):
        """Test that HTML plot files are created"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Check that HTML files were created
            html_files = list((output_dir / "raw").glob("*.html"))
            assert len(html_files) > 0, "No HTML files were created"


    def test_creates_single_html_with_all_plots(self, sample_results, sample_config):
        """Test that a single HTML file is created"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            html_files = list((output_dir / "raw").glob("*.html"))

            # Should create exactly one HTML file
            assert len(html_files) == 1, f"Expected 1 HTML file, found {len(html_files)}"

            # File should not be empty
            html_content = html_files[0].read_text()
            assert len(html_content) > 0, "HTML file is empty"

            # Should contain Plotly content
            assert 'plotly' in html_content.lower() or 'Plotly' in html_content

    def test_html_contains_multiple_plots(self, sample_results, sample_config):
        """Test that HTML contains multiple plot divs (one per chromosome)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            html_files = list((output_dir / "raw").glob("*.html"))
            html_content = html_files[0].read_text()

            # Count plot divs (Plotly creates divs with specific IDs)
            div_count = html_content.count('<div class=plot>')

            # Should have at least 2 plots (chr1 and chr2)
            assert div_count >= 2, f"Expected at least 2 plot divs, found {div_count}"

    def test_html_contains_all_mapping_modes(self, sample_results, sample_config):
        """Test that HTML contains plots for all mapping modes"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            # Add multiple mapping modes
            sample_results['fiveprime'] = copy.deepcopy(sample_results['threeprime'])

            create_metagene_plots(sample_results, sample_config, output_dir)

            html_files = list((output_dir / "raw").glob("*.html"))
            assert len(html_files) == 1

            html_content = html_files[0].read_text()

            # Should contain both mapping modes
            assert "3' Mapping" in html_content
            assert "5' Mapping" in html_content

    def test_html_file_naming(self, sample_results, sample_config):
        """Test that HTML file has expected naming pattern"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            create_metagene_plots(sample_results, sample_config, output_dir)

            html_files = list((output_dir / "raw").glob("*.html"))
            assert len(html_files) == 1

            # Check the filename (adjust based on your actual naming scheme)
            html_name = html_files[0].name

            assert html_name.endswith('.html')
            assert 'interactive_metagene_profiling' in html_name

    def test_creates_separate_html_per_normalization(self, sample_results, sample_config):
        """Test that separate HTML files are created for each normalization method"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"
            sample_config['normalizationMethods'] = ['raw', 'cpm']

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Should have one HTML in each normalization directory
            raw_html = list((output_dir / "raw").glob("*.html"))
            cpm_html = list((output_dir / "cpm").glob("*.html"))

            assert len(raw_html) == 1, f"Expected 1 HTML in raw/, found {len(raw_html)}"
            assert len(cpm_html) == 1, f"Expected 1 HTML in cpm/, found {len(cpm_html)}"


    def test_handles_multiple_normalization_methods(self, sample_results, sample_config):
        """Test that multiple normalization methods create separate directories"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"
            sample_config['normalizationMethods'] = ['raw', 'cpm']

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Check both directories exist
            assert (output_dir / "raw").exists()
            assert (output_dir / "cpm").exists()

            # Both should have Excel files
            assert (output_dir / "raw" / "threeprime_readcounts_start.xlsx").exists()
            assert (output_dir / "cpm" / "threeprime_readcounts_start.xlsx").exists()

    def test_handles_multiple_mapping_modes(self, sample_results, sample_config):
        """Test that multiple mapping modes are processed"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            # Add another mapping mode (deep copy of unique data)
            sample_results['fiveprime'] = copy.deepcopy(sample_results['threeprime'])

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Check Excel files for both modes exist
            assert (output_dir / "raw" / "threeprime_readcounts_start.xlsx").exists()
            assert (output_dir / "raw" / "fiveprime_readcounts_start.xlsx").exists()

    def test_cpm_normalization_applied(self, sample_results, sample_config):
        """Test that CPM normalization is correctly applied"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"
            sample_config['normalizationMethods'] = ['raw', 'cpm']

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Read both raw and CPM Excel files
            raw_path = output_dir / "raw" / "threeprime_readcounts_start.xlsx"
            cpm_path = output_dir / "cpm" / "threeprime_readcounts_start.xlsx"

            raw_df = pd.read_excel(raw_path, sheet_name='chr1')
            cpm_df = pd.read_excel(cpm_path, sheet_name='chr1')

            # CPM values should be different from raw
            # For read_length 28 with total_counts 10000:
            # raw value 150 -> CPM = (150 / 10000) * 1000000 = 15000
            pos_0_raw = raw_df[raw_df['position'] == 0]
            pos_0_cpm = cpm_df[cpm_df['position'] == 0]

            # Values should be different (normalized)
            if '28' in raw_df.columns:
                assert pos_0_raw['28'].iloc[0] != pos_0_cpm['28'].iloc[0]

    def test_handles_empty_results(self, sample_config):
        """Test handling of empty results dictionary"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"
            empty_results = {
                'threeprime': {
                    'start': {},
                    'stop': {},
                    'total_counts': {}
                }
            }

            # Should not crash
            create_metagene_plots(empty_results, sample_config, output_dir)

            # Directory should still be created
            assert output_dir.exists()

    def test_filters_read_lengths_for_plotting(self, sample_results, sample_config):
        """Test that only requested read lengths are plotted"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"
            # Only request read length 28
            sample_config['readLengths'] = '28-28'

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Excel should still contain all read lengths
            excel_path = output_dir / "raw" / "threeprime_readcounts_start.xlsx"
            df_dict = pd.read_excel(excel_path, sheet_name=None)
            chr1_df = df_dict['chr1']

            # Should have both 28 and 30 in Excel (all read lengths saved)
            assert '28' in chr1_df.columns
            assert '30' in chr1_df.columns

            html_files = list((output_dir / "raw").glob("*.html"))
            assert len(html_files) == 1

            html_content = html_files[0].read_text()

            assert '"name":"28 nt"' in html_content
            assert '"name":"30 nt"' not in html_content

    def test_multiple_output_formats(self, sample_results, sample_config):
        """Test that multiple output formats are created"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"
            sample_config['outputFormats'] = ['interactive', 'png']

            create_metagene_plots(sample_results, sample_config, output_dir)

            # Check for both HTML and PNG files
            html_files = list((output_dir / "raw").glob("*.html"))
            png_files = list((output_dir / "raw").glob("*.png"))

            assert len(html_files) > 0, "No HTML files created"
            # PNG creation depends on kaleido being installed
            assert len(png_files) > 0, "No PNG files created"


class TestExtractReadCountsIntegration:
    """Integration tests for extract_read_counts workflow"""

    @pytest.fixture
    def sample_annotation_file(self):
        """Create a minimal GFF/GTF annotation file for testing"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as tmp:
            # Write minimal GFF3 format annotation
            tmp.write("##gff-version 3\n")
            tmp.write("chr1\t.\tgene\t200\t499\t.\t+\t0\tID=gene1;Name=GenusOne\n")
            tmp.write("chr1\t.\tCDS\t200\t499\t.\t+\t0\tID=cds1;Parent=gene1\n")
            tmp.write("chr1\t.\tgene\t600\t899\t.\t-\t0\tID=gene2;Name=GenusTwo\n")
            tmp.write("chr1\t.\tCDS\t600\t899\t.\t+\t0\tID=cds2;Parent=gene2\n")
            tmp.write("chr1\t.\tgene\t1000\t1299\t.\t-\t0\tID=gene3;Name=GenusThree\n")
            tmp.write("chr1\t.\tCDS\t1000\t1299\t.\t-\t0\tID=cds3;Parent=gene3\n")
            tmp_path = tmp.name

        yield tmp_path

        # Cleanup
        Path(tmp_path).unlink(missing_ok=True)

    @pytest.fixture
    def sample_bam_file(self):
        """Create a minimal BAM file for testing"""
        with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as tmp:
            bam_path = tmp.name

        # Create BAM header
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 10000, 'SN': 'chr1'}]
        }

        # Create BAM file with some reads
        with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
            # Create aligned reads at various positions
            for i in range(50):
                a = pysam.AlignedSegment()
                a.query_name = f"read_{i}"
                a.query_sequence = "ATCG" * 7  # 28 nt
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 100 + (i * 5)  # Spread reads across positions
                a.mapping_quality = 60
                a.cigar = [(0, 28)]  # 28M
                a.query_qualities = pysam.qualitystring_to_array("I" * 28)
                a.set_tag('NH', 1)
                outf.write(a)

        # Index the BAM file
        pysam.index(bam_path)

        yield bam_path

        # Cleanup
        Path(bam_path).unlink(missing_ok=True)
        Path(f"{bam_path}.bai").unlink(missing_ok=True)

    @pytest.fixture
    def sample_config(self, sample_annotation_file):
        """Create sample configuration"""
        return {
            'annotationFilePath': sample_annotation_file,
            'mappingMethods': ['fiveprime'],
            'readLengths': '26-30',
            'positionsOutsideORF': 100,
            'positionsInsideORF': 250,
        }

    @pytest.fixture
    def genome_length_dict(self):
        """Sample genome lengths"""
        return {'chr1': 10000}

    def test_returns_results_dictionary(self, sample_bam_file, sample_config, genome_length_dict):
        """Test that function returns properly structured results"""
        results = extract_read_counts(
            sample_bam_file,
            "test_sample",
            sample_config,
            genome_length_dict
        )

        assert isinstance(results, dict)
        assert 'fiveprime' in results
        assert 'start' in results['fiveprime']
        assert 'stop' in results['fiveprime']
        assert 'total_counts' in results['fiveprime']

    def test_result_structure_matches_expected(self, sample_bam_file, sample_config, genome_length_dict):
        """Test that result structure matches expected nested dict format"""
        results = extract_read_counts(
            sample_bam_file,
            "test_sample",
            sample_config,
            genome_length_dict
        )

        start_dict = results['fiveprime']['start']
        stop_dict = results['fiveprime']['stop']
        total_counts = results['fiveprime']['total_counts']

        # All should be dictionaries
        assert isinstance(start_dict, dict)
        assert isinstance(stop_dict, dict)
        assert isinstance(total_counts, dict)

        # total_counts should have sample name
        assert 'test_sample' in total_counts
        assert isinstance(total_counts['test_sample'], dict)

    def test_processes_all_mapping_methods(self, sample_bam_file, sample_config, genome_length_dict):
        """Test that all specified mapping methods are processed"""
        sample_config['mappingMethods'] = ['fiveprime', 'threeprime']

        results = extract_read_counts(
            sample_bam_file,
            "test_sample",
            sample_config,
            genome_length_dict
        )

        assert 'fiveprime' in results
        assert 'threeprime' in results

        # Both should have the same structure
        for method in ['fiveprime', 'threeprime']:
            assert 'start' in results[method]
            assert 'stop' in results[method]
            assert 'total_counts' in results[method]

    def test_total_counts_includes_sample_name(self, sample_bam_file, sample_config, genome_length_dict):
        """Test that total_counts dictionary uses correct sample name"""
        sample_name = "my_custom_sample"

        results = extract_read_counts(
            sample_bam_file,
            sample_name,
            sample_config,
            genome_length_dict
        )

        total_counts = results['fiveprime']['total_counts']
        assert sample_name in total_counts

        # total_counts should map to read_length: count dict
        sample_counts = total_counts[sample_name]
        assert isinstance(sample_counts, dict)
        # Should have counts for different read lengths
        assert len(sample_counts) > 0

    def test_handles_empty_bam(self, sample_config, genome_length_dict):
        """Test handling of BAM file with no reads"""
        with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as tmp:
            bam_path = tmp.name

        try:
            # Create empty BAM
            header = {
                'HD': {'VN': '1.0'},
                'SQ': [{'LN': 10000, 'SN': 'chr1'}]
            }
            with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
                pass  # No reads

            pysam.index(bam_path)

            results = extract_read_counts(
                bam_path,
                "empty_sample",
                sample_config,
                genome_length_dict
            )

            # Should return structure even if empty/minimal
            assert isinstance(results, dict)
            assert 'fiveprime' in results

        finally:
            Path(bam_path).unlink(missing_ok=True)
            Path(f"{bam_path}.bai").unlink(missing_ok=True)

    def test_start_and_stop_dicts_exist(self, sample_bam_file, sample_config, genome_length_dict):
        """Test that start and stop dictionaries are created"""
        results = extract_read_counts(
            sample_bam_file,
            "test_sample",
            sample_config,
            genome_length_dict
        )

        start_dict = results['fiveprime']['start']
        stop_dict = results['fiveprime']['stop']

        # Both should be dictionaries
        assert isinstance(start_dict, dict)
        assert isinstance(stop_dict, dict)

    def test_prints_progress_information(self, sample_bam_file, sample_config, genome_length_dict, capsys):
        """Test that function prints progress information"""
        extract_read_counts(
            sample_bam_file,
            "test_sample",
            sample_config,
            genome_length_dict
        )

        captured = capsys.readouterr()

        # Should print sample name
        assert "test_sample" in captured.out

        # Should print mapping mode
        assert "fiveprime" in captured.out or "Mapping mode" in captured.out

        # Should print window counts
        assert "windows" in captured.out or "Found" in captured.out


class TestExtractReadCountsErrorHandling:
    """Tests for error handling in extract_read_counts"""

    @pytest.fixture
    def sample_annotation_file(self):
        """Create a minimal GFF/GTF annotation file for testing"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as tmp:
            # Write minimal GFF3 format annotation
            tmp.write("##gff-version 3\n")
            tmp.write("chr1\t.\tgene\t200\t499\t.\t+\t0\tID=gene1;Name=GenusOne\n")
            tmp.write("chr1\t.\tCDS\t200\t499\t.\t+\t0\tID=cds1;Parent=gene1\n")
            tmp.write("chr1\t.\tgene\t600\t899\t.\t-\t0\tID=gene2;Name=GenusTwo\n")
            tmp.write("chr1\t.\tCDS\t600\t899\t.\t+\t0\tID=cds2;Parent=gene2\n")
            tmp.write("chr1\t.\tgene\t1000\t1299\t.\t-\t0\tID=gene3;Name=GenusThree\n")
            tmp.write("chr1\t.\tCDS\t1000\t1299\t.\t-\t0\tID=cds3;Parent=gene3\n")
            tmp_path = tmp.name

        yield tmp_path

        # Cleanup
        Path(tmp_path).unlink(missing_ok=True)

    @pytest.fixture
    def sample_config(self, sample_annotation_file):
        """Create sample configuration"""
        return {
            'annotationFilePath': sample_annotation_file,
            'mappingMethods': ['fiveprime'],
            'readLengths': '26-30',
            'positionsOutsideORF': 100,
            'positionsInsideORF': 250,
        }


    @pytest.fixture
    def genome_length_dict(self):
        return {'chr1': 10000}

    def test_handles_missing_bam_file(self, sample_config, genome_length_dict):
        """Test handling of non-existent BAM file"""
        with pytest.raises(SystemExit, match="Error: Unable to read alignment file."):
            extract_read_counts(
                "nonexistent.bam",
                "test_sample",
                sample_config,
                genome_length_dict
            )

    def test_handles_invalid_bam_file(self, sample_config, genome_length_dict):
        """Test handling of invalid BAM file"""
        with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as tmp:
            tmp.write(b"not a valid bam file")
            tmp_path = tmp.name

        try:
            with pytest.raises(SystemExit, match="Error: Invalid BAM file format!"):  # pysam will raise various exceptions
                extract_read_counts(
                    tmp_path,
                    "test_sample",
                    sample_config,
                    genome_length_dict
                )
        finally:
            Path(tmp_path).unlink(missing_ok=True)

    def test_handles_missing_bam_index(self, sample_config, genome_length_dict):
        """Test handling of BAM file without index"""
        with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as tmp:
            bam_path = tmp.name

        try:
            # Create a valid BAM file with multiple reads
            # (more reads = more likely to need index)
            header = {
                'HD': {'VN': '1.0'},
                'SQ': [{'LN': 10000, 'SN': 'chr1'}]
            }
            with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
                for i in range(100):  # Create many reads
                    a = pysam.AlignedSegment()
                    a.query_name = f"read_{i}"
                    a.query_sequence = "ATCG" * 7
                    a.flag = 0
                    a.reference_id = 0
                    a.reference_start = 100 + (i * 10)
                    a.mapping_quality = 60
                    a.cigar = [(0, 28)]
                    a.query_qualities = pysam.qualitystring_to_array("I" * 28)
                    a.set_tag('NH', 1)
                    outf.write(a)

            # Remove index if it was created
            for index_path in [f"{bam_path}.bai", f"{bam_path[:-4]}.bai"]:
                if Path(index_path).exists():
                    Path(index_path).unlink()

            # Should fail with index error
            with pytest.raises(SystemExit, match="appropriate index file|index"):
                extract_read_counts(
                    bam_path,
                    "test_sample",
                    sample_config,
                    genome_length_dict
                )

        finally:
            Path(bam_path).unlink(missing_ok=True)
            for index_path in [f"{bam_path}.bai", f"{bam_path[:-4]}.bai"]:
                Path(index_path).unlink(missing_ok=True)