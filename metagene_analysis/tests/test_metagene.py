import pytest
from unittest.mock import patch
import interlap


from lib.metagene import (
    extract_start_metagene_windows,
    extract_stop_metagene_windows,
    create_interlap_dict
)

class TestExtractStartMetageneWindows:
    """Test suite for extract_start_metagene_windows function"""

    @pytest.fixture
    def basic_config(self):
        """Basic configuration"""
        return {
            "positionsOutsideORF": 50,
            "positionsInsideORF": 100
        }

    @pytest.fixture
    def annotation_identifiers(self):
        """Sample annotation identifiers (chrom, start, stop, strand)"""
        return [
            ("chr1", 1000, 1300, "+"),
            ("chr1", 2000, 2500, "-"),
            ("chr2", 500, 800, "+")
        ]

    def test_returns_list(self, annotation_identifiers, basic_config):
        """Test that function returns a list"""
        result = extract_start_metagene_windows(annotation_identifiers, basic_config)
        assert isinstance(result, list)

    def test_window_count(self, annotation_identifiers, basic_config):
        """Test that function returns correct number of windows"""
        result = extract_start_metagene_windows(annotation_identifiers, basic_config)
        assert len(result) == 3

    def test_positive_strand_coordinates(self, basic_config):
        """Test coordinate calculation for positive strand"""
        identifiers = [("chr1", 1000, 1300, "+")]
        result = extract_start_metagene_windows(identifiers, basic_config)

        _, upstream, downstream, strand, _ = result[0]

        # upstream_region = start + window[0] = 1000 + (-50) = 950
        # downstream_region = start + window[1] = 1000 + 100 = 1100
        assert upstream == 950
        assert downstream == 1100
        assert strand == "+"

    def test_negative_strand_coordinates(self, basic_config):
        """Test coordinate calculation for negative strand"""
        identifiers = [("chr1", 2000, 2500, "-")]
        result = extract_start_metagene_windows(identifiers, basic_config)

        _, upstream, downstream, strand, _ = result[0]

        # upstream_region = stop - window[1] = 2500 - 100 = 2400
        # downstream_region = stop - window[0] = 2500 - (-50) = 2550
        assert upstream == 2400
        assert downstream == 2550
        assert strand == "-"

    def test_output_id_format(self, basic_config):
        """Test that output ID is formatted correctly (1-based)"""
        identifiers = [("chr1", 1000, 1300, "+")]
        result = extract_start_metagene_windows(identifiers, basic_config)

        _, _, _, _, out_id = result[0]
        assert out_id == "chr1:1001-1301:+"

    def test_skips_cds_shorter_than_window(self):
        """Test that short CDSs are skipped when window extends beyond CDS"""
        config = {
            "positionsOutsideORF": 50,
            "positionsInsideORF": 100
        }
        # CDS length = 50, but positionsInsideORF = 100 (longer than CDS)
        identifiers = [("chr1", 1000, 1050, "+")]
        result = extract_start_metagene_windows(identifiers, config)

        assert len(result) == 0

    def test_includes_cds_longer_than_window(self):
        """Test that CDSs longer than window are included"""
        config = {
            "positionsOutsideORF": 50,
            "positionsInsideORF": 100
        }
        # CDS length = 300, longer than window
        identifiers = [("chr1", 1000, 1300, "+")]
        result = extract_start_metagene_windows(identifiers, config)

        assert len(result) == 1

    def test_window_validation_error(self):
        """Test that invalid window configuration causes exit"""
        config = {
            "positionsOutsideORF": -100,  # invalid negative value
            "positionsInsideORF": 50
        }
        identifiers = [("chr1", 1000, 1300, "+")]

        with pytest.raises(SystemExit):
            extract_start_metagene_windows(identifiers, config)

    def test_empty_input(self, basic_config):
        """Test with empty annotation list"""
        result = extract_start_metagene_windows([], basic_config)
        assert not result

    def test_tuple_structure(self, annotation_identifiers, basic_config):
        """Test that each window is a tuple with correct structure"""
        result = extract_start_metagene_windows(annotation_identifiers, basic_config)

        for window in result:
            assert isinstance(window, tuple)
            assert len(window) == 5
            chrom, upstream, downstream, strand, out_id = window
            assert isinstance(chrom, str)
            assert isinstance(upstream, int)
            assert isinstance(downstream, int)
            assert strand in ["+", "-"]
            assert isinstance(out_id, str)


class TestExtractStopMetageneWindows:
    """Test suite for extract_stop_metagene_windows function"""

    @pytest.fixture
    def basic_config(self):
        """Basic configuration"""
        return {
            "positionsOutsideORF": 50,
            "positionsInsideORF": 100
        }

    @pytest.fixture
    def annotation_identifiers(self):
        """Sample annotation identifiers"""
        return [
            ("chr1", 1000, 1300, "+"),
            ("chr1", 2000, 2500, "-"),
        ]

    def test_returns_list(self, annotation_identifiers, basic_config):
        """Test that function returns a list"""
        result = extract_stop_metagene_windows(annotation_identifiers, basic_config)
        assert isinstance(result, list)

    def test_positive_strand_coordinates(self, basic_config):
        """Test coordinate calculation for positive strand"""
        identifiers = [("chr1", 1000, 1300, "+")]
        result = extract_stop_metagene_windows(identifiers, basic_config)

        chrom, upstream, downstream, strand, out_id = result[0]

        # window = [-100, 50]
        # upstream_region = stop + window[0] = 1300 + (-100) = 1200
        # downstream_region = stop + window[1] = 1300 + 50 = 1350
        assert upstream == 1200
        assert downstream == 1350
        assert strand == "+"

    def test_negative_strand_coordinates(self, basic_config):
        """Test coordinate calculation for negative strand"""
        identifiers = [("chr1", 2000, 2500, "-")]
        result = extract_stop_metagene_windows(identifiers, basic_config)

        chrom, upstream, downstream, strand, out_id = result[0]

        # window = [-100, 50]
        # upstream_region = start - window[1] = 2000 - 50 = 1950
        # downstream_region = start - window[0] = 2000 - (-100) = 2100
        assert upstream == 1950
        assert downstream == 2100
        assert strand == "-"

    def test_skips_cds_shorter_than_inside_window(self):
        """Test that short CDSs are skipped when inside window extends beyond CDS"""
        config = {
            "positionsOutsideORF": 50,
            "positionsInsideORF": 100
        }
        # CDS length = 50, but positionsInsideORF = 100
        identifiers = [("chr1", 1000, 1050, "+")]
        result = extract_stop_metagene_windows(identifiers, config)

        assert len(result) == 0

    def test_includes_cds_longer_than_window(self):
        """Test that CDSs longer than window are included"""
        config = {
            "positionsOutsideORF": 50,
            "positionsInsideORF": 100
        }
        identifiers = [("chr1", 1000, 1300, "+")]
        result = extract_stop_metagene_windows(identifiers, config)

        assert len(result) == 1

    def test_window_validation_error(self):
        """Test that invalid window configuration causes exit"""
        config = {
            "positionsOutsideORF": 50,
            "positionsInsideORF": -100
        }
        identifiers = [("chr1", 1000, 1300, "+")]

        with pytest.raises(SystemExit):
            extract_stop_metagene_windows(identifiers, config)

    def test_empty_input(self, basic_config):
        """Test with empty annotation list"""
        result = extract_stop_metagene_windows([], basic_config)
        assert not result


class TestCreateInterlapDict:
    """Test suite for create_interlap_dict function"""

    @pytest.fixture
    def basic_windows(self):
        """Sample metagene windows"""
        return [
            ("chr1", 100, 200, "+", "gene1"),
            ("chr1", 300, 400, "+", "gene2"),
            ("chr1", 500, 600, "-", "gene3"),
            ("chr2", 100, 200, "+", "gene4")
        ]

    @pytest.fixture
    def empty_locus_map(self):
        """Empty locus map"""
        return {}

    def test_returns_dict(self, basic_windows, empty_locus_map):
        """Test that function returns a dictionary"""
        result = create_interlap_dict(basic_windows, empty_locus_map)
        assert isinstance(result, dict)

    def test_keys_are_chrom_strand_tuples(self, basic_windows, empty_locus_map):
        """Test that keys are (chromosome, strand) tuples"""
        result = create_interlap_dict(basic_windows, empty_locus_map)

        for key in result:
            assert isinstance(key, tuple)
            assert len(key) == 2
            chrom, strand = key
            assert isinstance(chrom, str)
            assert strand in ["+", "-"]

    def test_correct_number_of_keys(self, basic_windows, empty_locus_map):
        """Test that we have correct number of (chrom, strand) combinations"""
        result = create_interlap_dict(basic_windows, empty_locus_map)

        # chr1 +, chr1 -, chr2 +
        assert len(result) == 3
        assert ("chr1", "+") in result
        assert ("chr1", "-") in result
        assert ("chr2", "+") in result
        assert ("chr2", "-") not in result

    def test_values_are_interlap_objects(self, basic_windows, empty_locus_map):
        """Test that values are InterLap objects"""
        result = create_interlap_dict(basic_windows, empty_locus_map)

        for value in result.values():
            assert isinstance(value, interlap.InterLap)

    def test_interlap_contains_correct_ranges(self, basic_windows, empty_locus_map):
        """Test that InterLap objects contain correct ranges"""
        result = create_interlap_dict(basic_windows, empty_locus_map)

        # Check chr1 + strand
        chr1_plus = result[("chr1", "+")]
        matches = list(chr1_plus.find((100, 200)))
        assert len(matches) == 1
        assert matches[0][2] == "gene1"

        matches = list(chr1_plus.find((300, 400)))
        assert len(matches) == 1
        assert matches[0][2] == "gene2"

    def test_locus_map_remapping(self):
        """Test that identifiers are remapped using locus_map"""
        windows = [
            ("chr1", 100, 200, "+", "gene1"),
            ("chr1", 300, 400, "+", "gene2")
        ]
        locus_map = {
            "gene1": "locus_A",
            "gene2": "locus_B"
        }

        result = create_interlap_dict(windows, locus_map)

        chr1_plus = result[("chr1", "+")]
        matches = list(chr1_plus.find((100, 200)))
        assert matches[0][2] == "locus_A"

        matches = list(chr1_plus.find((300, 400)))
        assert matches[0][2] == "locus_B"

    def test_locus_map_partial_remapping(self):
        """Test that only mapped identifiers are remapped"""
        windows = [
            ("chr1", 100, 200, "+", "gene1"),
            ("chr1", 300, 400, "+", "gene2")
        ]
        locus_map = {
            "gene1": "locus_A",
            # gene2 not in map
        }

        result = create_interlap_dict(windows, locus_map)

        chr1_plus = result[("chr1", "+")]
        matches = list(chr1_plus.find((100, 200)))
        assert matches[0][2] == "locus_A"

        matches = list(chr1_plus.find((300, 400)))
        assert matches[0][2] == "gene2"  # unchanged

    def test_empty_windows(self, empty_locus_map):
        """Test with empty windows list"""
        result = create_interlap_dict([], empty_locus_map)
        assert not result

    def test_overlapping_windows_same_strand(self, empty_locus_map):
        """Test handling of overlapping windows on same strand"""
        windows = [
            ("chr1", 100, 250, "+", "gene1"),
            ("chr1", 200, 300, "+", "gene2")
        ]

        result = create_interlap_dict(windows, empty_locus_map)

        chr1_plus = result[("chr1", "+")]
        # Both should be findable in the overlap region
        matches = list(chr1_plus.find((200, 250)))
        assert len(matches) == 2

    def test_same_coordinates_different_strands(self, empty_locus_map):
        """Test that same coordinates on different strands are kept separate"""
        windows = [
            ("chr1", 100, 200, "+", "gene1"),
            ("chr1", 100, 200, "-", "gene2")
        ]

        result = create_interlap_dict(windows, empty_locus_map)

        assert ("chr1", "+") in result
        assert ("chr1", "-") in result

        matches_plus = list(result[("chr1", "+")].find((100, 200)))
        matches_minus = list(result[("chr1", "-")].find((100, 200)))

        assert len(matches_plus) == 1
        assert len(matches_minus) == 1
        assert matches_plus[0][2] == "gene1"
        assert matches_minus[0][2] == "gene2"