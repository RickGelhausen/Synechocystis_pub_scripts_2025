"""
Unit tests for lib.annotation module

Tests for:
- create_annotation_intervals_dict
- filter_annotation_identifiers
"""

import pytest
import pandas as pd
import interlap
from unittest.mock import patch
from io import StringIO

from lib.annotation import create_annotation_intervals_dict
from lib.annotation import filter_annotation_identifiers

class TestCreateAnnotationIntervalsDict:
    """Test suite for create_annotation_intervals_dict function"""

    @pytest.fixture
    def simple_annotation_df(self):
        """Create a simple annotation DataFrame for testing"""
        data = {
            0: ["chr1", "chr1", "chr1", "chr2"],  # chromosome
            1: ["source", "source", "source", "source"],  # source
            2: ["CDS", "CDS", "gene", "CDS"],  # feature type
            3: [101, 201, 301, 401],  # start (1-based in GFF)
            4: [200, 300, 400, 500],  # end (1-based in GFF)
            5: [".", ".", ".", "."],  # score
            6: ["+", "+", "+", "-"],  # strand
            7: [".", ".", ".", "."],  # frame
            8: ["ID=cds1", "ID=cds2", "ID=gene1", "ID=cds3"],  # attributes
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def mixed_features_df(self):
        """DataFrame with CDS and non-CDS features"""
        data = {
            0: ["chr1", "chr1", "chr1", "chr1"],
            1: ["source"] * 4,
            2: ["CDS", "gene", "exon", "CDS"],  # Mix of features
            3: [101, 201, 301, 401],
            4: [200, 300, 400, 500],
            5: ["."] * 4,
            6: ["+", "+", "+", "+"],
            7: ["."] * 4,
            8: ["ID=cds1", "ID=gene1", "ID=exon1", "ID=cds2"],
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def empty_annotation_df(self):
        """Empty annotation DataFrame"""
        columns = list(range(9))
        return pd.DataFrame(columns=columns)

    def test_creates_dict_with_correct_structure(self, simple_annotation_df):
        """Test that function returns a dictionary with correct structure"""


        result = create_annotation_intervals_dict(simple_annotation_df)

        assert isinstance(result, dict)
        assert len(result) > 0

        # Check that keys are tuples of (chromosome, strand)
        for key in result:
            assert isinstance(key, tuple)
            assert len(key) == 2
            assert isinstance(key[0], str)  # chromosome
            assert key[1] in ["+", "-"]  # strand

    def test_creates_interlap_objects(self, simple_annotation_df):
        """Test that dictionary values are InterLap objects"""

        result = create_annotation_intervals_dict(simple_annotation_df)

        for value in result.values():
            assert isinstance(value, interlap.InterLap)

    def test_filters_only_cds_features(self, mixed_features_df):
        """Test that only CDS features are included"""

        result = create_annotation_intervals_dict(mixed_features_df)

        # Should only have chr1, + strand (2 CDS features there)
        assert ("chr1", "+") in result

        # Get intervals and check count
        intervals = list(result[("chr1", "+")].find((0, 1000)))
        assert len(intervals) == 2  # Only 2 CDS features
        assert intervals == [(100, 199), (400, 499)]  # Correct intervals

    def test_converts_coordinates_to_zero_based(self, simple_annotation_df):
        """Test that 1-based GFF coordinates are converted to 0-based"""

        result = create_annotation_intervals_dict(simple_annotation_df)

        # Get intervals for chr1, plus strand
        intervals = list(result[("chr1", "+")].find((0, 1000)))
        intervals.sort(key=lambda x: x[0])

        # GFF has start=101, end=200 (1-based, inclusive)
        # Should become 100, 199 (0-based, exclusive end becomes inclusive)
        assert intervals[0] == (100, 199)  # First CDS
        assert intervals[1] == (200, 299)  # Second CDS

    def test_groups_by_chromosome_and_strand(self, simple_annotation_df):
        """Test that features are correctly grouped by chromosome and strand"""

        result = create_annotation_intervals_dict(simple_annotation_df)

        # Should have entries for chr1+ and chr2-
        assert ("chr1", "+") in result
        assert ("chr2", "-") in result

        # Check counts
        intervals_chr1_plus = list(result[("chr1", "+")].find((0, 1000)))
        intervals_chr2_minus = list(result[("chr2", "-")].find((0, 1000)))

        assert len(intervals_chr1_plus) == 2  # Two CDS on chr1+
        assert len(intervals_chr2_minus) == 1  # One CDS on chr2-

    def test_handles_empty_dataframe(self, empty_annotation_df):
        """Test that function handles empty DataFrame gracefully"""

        result = create_annotation_intervals_dict(empty_annotation_df)

        assert isinstance(result, dict)
        assert len(result) == 0

    def test_case_insensitive_cds_filtering(self):
        """Test that CDS filtering is case-insensitive"""

        # Create DataFrame with different case variations
        data = {
            0: ["chr1", "chr1", "chr1"],
            1: ["source"] * 3,
            2: ["CDS", "cds", "Cds"],  # Different cases
            3: [101, 201, 301],
            4: [200, 300, 400],
            5: ["."] * 3,
            6: ["+"] * 3,
            7: ["."] * 3,
            8: ["ID=1", "ID=2", "ID=3"],
        }
        df = pd.DataFrame(data)

        result = create_annotation_intervals_dict(df)

        intervals = list(result[("chr1", "+")].find((0, 1000)))
        assert len(intervals) == 3  # All three CDS features included

    def test_intervals_can_be_queried(self, simple_annotation_df):
        """Test that created InterLap objects can be queried"""

        result = create_annotation_intervals_dict(simple_annotation_df)

        # Query for overlaps at position 150 (should hit first CDS: 100-199)
        overlaps = list(result[("chr1", "+")].find((150, 150)))
        assert len(overlaps) == 1
        assert overlaps[0] == (100, 199)

        # Query for overlaps at position 250 (should hit second CDS: 200-299)
        overlaps = list(result[("chr1", "+")].find((250, 250)))
        assert len(overlaps) == 1
        assert overlaps[0] == (200, 299)

        # Query gap (should hit nothing)
        overlaps = list(result[("chr1", "+")].find((50, 50)))
        assert len(overlaps) == 0


class TestFilterAnnotationIdentifiers:
    """Test suite for filter_annotation_identifiers function"""

    @pytest.fixture
    def basic_config(self):
        """Basic configuration for testing"""
        return {
            "annotationFilePath": "test_annotation.gff",
            "filteringMethods": [],
            "rpkmThreshold": 5.0,
            "neighboringGenesDistance": 100,
            "lengthCutoff": 50,
            "positionsOutsideORF": 50,
            "positionsInsideORF": 100,
        }

    @pytest.fixture
    def mock_annotation_content(self):
        """Mock annotation file content"""
        return """chr1\tsource\tgene\t1\t300\t.\t+\t.\tID=gene1
chr1\tsource\tCDS\t101\t301\t.\t+\t.\tParent=gene1;locus_tag=TAG001
chr1\tsource\tCDS\t401\t601\t.\t+\t.\tlocus_tag=TAG002
chr2\tsource\tpseudogene\t1\t300\t.\t+\t.\tID=gene2
chr2\tsource\tCDS\t101\t301\t.\t+\t.\tParent=gene2;locus_tag=TAG003
"""

    @pytest.fixture
    def mock_read_intervals_dict(self):
        """Mock read intervals dictionary"""
        intervals_dict = {}
        inter = interlap.InterLap()
        inter.update([(100, 299)])  # Reads covering the gene
        intervals_dict[("chr1", "+")] = inter
        return intervals_dict

    @pytest.fixture
    def mock_total_counts_dict(self):
        """Mock total read counts per chromosome"""
        return {"chr1": 1000000, "chr2": 500000}

    @pytest.fixture
    def mock_genome_length_dict(self):
        """Mock genome lengths"""
        return {"chr1": 10000, "chr2": 10000}

    def test_returns_correct_structure(self, basic_config, mock_annotation_content):
        """Test that function returns tuple of (identifiers, locus_map)"""

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(mock_annotation_content), sep="\t", comment="#", header=None
        )

        # Patch where pd.read_csv is USED (in lib.annotation), not where it's defined
        with patch('lib.annotation.pd.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, locus_map = filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

        assert isinstance(identifiers, list)
        assert isinstance(locus_map, dict)


    def test_filters_pseudogenes(self, basic_config, mock_annotation_content):
        """Test that pseudogenes are filtered out"""

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(mock_annotation_content), sep="\t", comment="#", header=None
        )

        with patch('lib.annotation.pd.read_csv', return_value=df) as mock_read:
            # Remove the print patch temporarily to see what the function prints
            # with patch('builtins.print'):
            identifiers, locus_map = filter_annotation_identifiers(
                None, None, None, "fiveprime", basic_config
            )

        assert len(identifiers) == 2
        assert "TAG001" in locus_map.values()
        assert "TAG002" in locus_map.values()
        assert "TAG003" not in locus_map.values()

    def test_filters_by_length(self, basic_config):
        """Test length-based filtering"""

        annotation_content = """chr1\tsource\tCDS\t101\t150\t.\t+\t.\tlocus_tag=TAG001
chr1\tsource\tCDS\t201\t500\t.\t+\t.\tlocus_tag=TAG002
"""
        basic_config["filteringMethods"] = ["length"]
        basic_config["lengthCutoff"] = 100

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="	", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, locus_map = filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

        # Only TAG002 should pass (length = 300)
        assert len(identifiers) == 1
        assert "TAG002" in locus_map.values()
        assert "TAG001" not in locus_map.values()

    def test_filters_by_overlap(self, basic_config):
        """Test overlap-based filtering"""

        # Two overlapping genes
        annotation_content = """chr1\tsource\tCDS\t101\t300\t.\t+\t.\tlocus_tag=TAG001
chr1\tsource\tCDS\t250\t500\t.\t+\t.\tlocus_tag=TAG002
"""
        basic_config["filteringMethods"] = ["overlap"]
        basic_config["neighboringGenesDistance"] = 0

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="	", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, _ = filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

        # Both genes overlap, so both should be filtered
        assert len(identifiers) == 0

    def test_filters_by_rpkm(self, basic_config, mock_read_intervals_dict,
                             mock_total_counts_dict):
        """Test RPKM-based filtering"""

        annotation_content = """chr1\tsource\tCDS\t101\t301\t.\t+\t.\tlocus_tag=TAG001
chr1\tsource\tCDS\t401\t601\t.\t+\t.\tlocus_tag=TAG002
"""
        basic_config["filteringMethods"] = ["rpkm"]
        basic_config["rpkmThreshold"] = 10.0

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="\t", comment="#", header=None
        )

        with patch('lib.annotation.pd.read_csv', return_value=df):
            with patch('lib.misc.count_reads') as mock_count:
                with patch('lib.misc.calculate_rpkm') as mock_rpkm:
                    # First gene has high RPKM, second has low
                    mock_count.side_effect = [100, 1]
                    mock_rpkm.side_effect = [50.0, 0.5]

                    with patch('builtins.print'):
                        identifiers, locus_map = filter_annotation_identifiers(
                            mock_read_intervals_dict,
                            mock_total_counts_dict,
                            None,
                            "fiveprime",
                            basic_config
                        )

        # Only TAG001 should pass (RPKM = 50.0 > threshold)
        assert len(identifiers) == 1
        assert "TAG001" in locus_map.values()

    def test_filters_by_boundary(self, basic_config, mock_genome_length_dict):
        """Test boundary-based filtering"""

        # Gene too close to chromosome start
        annotation_content = """chr1\tsource\tCDS\t1\t201\t.\t+\t.\tlocus_tag=TAG001
chr1\tsource\tCDS\t501\t701\t.\t+\t.\tlocus_tag=TAG002
"""
        basic_config["positionsOutsideORF"] = 100

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="\t", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, locus_map = filter_annotation_identifiers(
                    None, None, mock_genome_length_dict, "fiveprime", basic_config
                )

        # TAG001 should be filtered (start - 100 < 0)
        # TAG002 should pass
        assert len(identifiers) == 1
        assert "TAG002" in locus_map.values()
        assert "TAG001" not in locus_map.values()

    def test_filters_non_divisible_by_three(self, basic_config):
        """Test that genes not divisible by 3 are filtered"""

        # First gene: length = 200 (divisible by 3? No: 200 % 3 = 2)
        # Second gene: length = 201 (divisible by 3? Yes: 201 % 3 = 0)
        annotation_content = """chr1\tsource\tCDS\t101\t300\t.\t+\t.\tlocus_tag=TAG001
chr1\tsource\tCDS\t401\t601\t.\t+\t.\tlocus_tag=TAG002
"""

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="\t", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, locus_map = filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

        # Only TAG002 should pass (length divisible by 3)
        assert len(identifiers) == 1
        assert "TAG002" in locus_map.values()

    def test_returns_zero_based_coordinates(self, basic_config):
        """Test that returned coordinates are 0-based"""

        annotation_content = """chr1\tsource\tCDS\t101\t301\t.\t+\t.\tlocus_tag=TAG001"""

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="	", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, _ = filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

        # GFF: 101-301 (1-based) -> should return 100-300 (0-based)
        assert identifiers[0] == ("chr1", 100, 300, "+")

    def test_locus_map_format(self, basic_config):
        """Test that locus_map has correct format"""

        annotation_content = """chr1\tsource\tCDS\t101\t301\t.\t+\t.\tlocus_tag=TAG001"""

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="\t", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                _, locus_map = filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

        # Locus map key format: "chrom:start-end:strand" (1-based for output)
        expected_key = "chr1:101-301:+"
        assert expected_key in locus_map
        assert locus_map[expected_key] == "TAG001"

    def test_multiple_filtering_methods(self, basic_config):
        """Test that multiple filtering methods work together"""

        annotation_content = """chr1\tsource\tCDS\t101\t201\t.\t+\t.\tlocus_tag=TAG001
chr1\tsource\tCDS\t301\t501\t.\t+\t.\tlocus_tag=TAG002
chr1\tsource\tCDS\t601\t801\t.\t+\t.\tlocus_tag=TAG003
chr1\tsource\tCDS\t810\t1001\t.\t+\t.\tlocus_tag=TAG004
"""
        basic_config["filteringMethods"] = ["length", "overlap"]
        basic_config["lengthCutoff"] = 150
        basic_config["neighboringGenesDistance"] = 50

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="\t", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, locus_map = filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

        # TAG001: length=100, filtered by length
        # TAG002: length=200, passes length filter
        # TAG003: length=200, passes length filter
        assert len(identifiers) == 1
        assert "TAG002" in locus_map.values()
        assert "TAG003" not in locus_map.values()

    def test_rpkm_filtering_raises_error_without_required_data(self, basic_config):
        """Test that RPKM filtering raises error when data is missing"""

        annotation_content = """chr1\tsource\tCDS\t101\t301\t.\t+\t.\tlocus_tag=TAG001"""
        basic_config["filteringMethods"] = ["rpkm"]

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="\t", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                with pytest.raises(ValueError, match="RPKM filtering requires"):
                    filter_annotation_identifiers(
                        None, None, None, "fiveprime", basic_config
                    )

    def test_prints_filtering_statistics(self, basic_config, mock_annotation_content):
        """Test that filtering statistics are printed"""

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(mock_annotation_content), sep="\t", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print') as mock_print:
                filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

                # Check that print was called with statistics
                assert mock_print.called
                # Verify some expected output
                calls = [str(call) for call in mock_print.call_args_list]
                stats_output = ''.join(calls)
                assert "Included genes" in stats_output
                assert "Excluded genes" in stats_output

    def test_handles_missing_locus_tag(self, basic_config):
        """Test that genes without locus_tag are still processed"""

        annotation_content = """chr1\tsource\tCDS\t101\t301\t.\t+\t.\tID=cds001
chr1\tsource\tCDS\t401\t601\t.\t+\t.\tlocus_tag=TAG002
"""

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="\t", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, locus_map = filter_annotation_identifiers(
                    None, None, None, "fiveprime", basic_config
                )

        # Both should be in identifiers
        assert len(identifiers) == 2
        # Only one should be in locus_map
        assert len(locus_map) == 1
        assert "TAG002" in locus_map.values()

    def test_boundary_filtering_minus_strand(self, basic_config,
                                             mock_genome_length_dict):
        """Test boundary filtering for minus strand genes"""

        # Gene on minus strand close to chromosome end
        annotation_content = """chr1\tsource\tCDS\t9902\t10000\t.\t-\t.\tlocus_tag=TAG001
chr1\tsource\tCDS\t5000\t5200\t.\t-\t.\tlocus_tag=TAG002
"""
        basic_config["positionsOutsideORF"] = 100
        mock_genome_length_dict["chr1"] = 10000

        # Create DataFrame BEFORE patching
        df = pd.read_csv(
            StringIO(annotation_content), sep="\t", comment="#", header=None
        )

        with patch('pandas.read_csv', return_value=df):
            with patch('builtins.print'):
                identifiers, locus_map = filter_annotation_identifiers(
                    None, None, mock_genome_length_dict, "fiveprime", basic_config
                )

        # TAG001 should be filtered (end + 100 would go past chromosome end)
        # TAG002 should pass
        assert len(identifiers) == 1
        assert "TAG002" in locus_map.values()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])