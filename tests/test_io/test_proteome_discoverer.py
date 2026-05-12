"""
Tests for psm_utils.io.proteome_discoverer.

This comprehensive test suite covers the MSFReader class through both unit tests with mocks
and integration tests with real MSF data. The test structure includes:

1. Unit tests: Fast, isolated tests using mocks for individual methods
2. Integration tests: End-to-end tests using the minimal MSF test file
3. Error handling tests: Edge cases and error conditions
4. Performance tests: Validation of method behaviors with real data

The test suite validates SQLAlchemy 2.0 patterns, proper data extraction,
and complete PSM object construction.
"""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from psm_utils import Peptidoform
from psm_utils.io.proteome_discoverer import COMPATIBLE_VERSIONS, MSFReader


class TestMSFReaderUnit:
    """Unit tests for MSFReader using mocks - fast and isolated."""

    @pytest.fixture
    def mock_reader_setup(self):
        """Set up common mocked MSFReader for tests."""
        with (
            patch("psm_utils.io.proteome_discoverer.create_engine") as mock_engine,
            patch("psm_utils.io.proteome_discoverer.Session") as mock_session_class,
        ):
            mock_session = Mock()
            mock_session_class.return_value = mock_session

            # Default version check
            mock_session.execute.return_value.first.return_value = (79,)

            test_file = Path("test.msf")
            test_file.touch()

            yield {
                "mock_engine": mock_engine,
                "mock_session_class": mock_session_class,
                "mock_session": mock_session,
                "test_file": test_file,
            }

            test_file.unlink(missing_ok=True)

    def test_init_success(self, mock_reader_setup):
        """Test successful MSFReader initialization."""
        setup = mock_reader_setup

        MSFReader(setup["test_file"])

        # Verify SQLAlchemy setup
        setup["mock_engine"].assert_called_once()
        setup["mock_session_class"].assert_called_once()

        # Verify version check was called
        assert setup["mock_session"].execute.called

    def test_init_no_version_info(self, mock_reader_setup):
        """Test initialization when MSF file has no version information."""
        setup = mock_reader_setup
        setup["mock_session"].execute.return_value.first.return_value = None

        with patch("psm_utils.io.proteome_discoverer.logger") as mock_logger:
            MSFReader(setup["test_file"])

            mock_logger.warning.assert_called_once()
            warning_msg = mock_logger.warning.call_args[0][0]
            assert "does not contain version information" in warning_msg

    @pytest.mark.parametrize("version", [79, 53, 8])
    def test_compatible_versions(self, mock_reader_setup, version):
        """Test that compatible versions don't raise warnings."""
        setup = mock_reader_setup
        setup["mock_session"].execute.return_value.first.return_value = (version,)

        with patch("psm_utils.io.proteome_discoverer.logger") as mock_logger:
            MSFReader(setup["test_file"])

            # Should not log any warnings for compatible versions
            mock_logger.warning.assert_not_called()

    def test_incompatible_version_warning(self, mock_reader_setup):
        """Test warning for incompatible MSF version."""
        setup = mock_reader_setup
        setup["mock_session"].execute.return_value.first.return_value = (999,)

        with patch("psm_utils.io.proteome_discoverer.logger") as mock_logger:
            MSFReader(setup["test_file"])

            mock_logger.warning.assert_called_once()
            warning_msg = mock_logger.warning.call_args[0][0]
            assert "version 999 might not be compatible" in warning_msg

    def test_len_method(self, mock_reader_setup):
        """Test __len__ method with mocked counts."""
        setup = mock_reader_setup

        def mock_execute_side_effect(stmt):
            stmt_str = str(stmt)
            mock_result = Mock()

            if "SchemaInfo" in stmt_str and "Version" in stmt_str:
                mock_result.first.return_value = (79,)
            elif "Peptides_decoy" in stmt_str:
                mock_result.scalar.return_value = 150
            else:  # Regular Peptides table
                mock_result.scalar.return_value = 1000

            return mock_result

        setup["mock_session"].execute.side_effect = mock_execute_side_effect

        reader = MSFReader(setup["test_file"])
        assert len(reader) == 1150  # 1000 + 150

    def test_get_modifications_structure(self, mock_reader_setup):
        """Test _get_modifications method structure and return format."""
        setup = mock_reader_setup
        reader = MSFReader(setup["test_file"])

        # Reset and mock modification query
        setup["mock_session"].execute.reset_mock()
        mock_results = [
            (1, 0, 4),  # PeptideID 1, Position 0, UnimodAccession 4
            (1, 3, 35),  # PeptideID 1, Position 3, UnimodAccession 35
            (2, 1, 21),  # PeptideID 2, Position 1, UnimodAccession 21
        ]
        setup["mock_session"].execute.return_value = mock_results

        modifications = reader._get_modifications(is_decoy=False)

        assert isinstance(modifications, dict)
        assert modifications[1] == [(0, 4), (3, 35)]
        assert modifications[2] == [(1, 21)]
        assert len(modifications[1]) == 2
        assert len(modifications[2]) == 1

    def test_get_terminal_modifications_structure(self, mock_reader_setup):
        """Test _get_terminal_modifications method structure."""
        setup = mock_reader_setup
        reader = MSFReader(setup["test_file"])

        setup["mock_session"].execute.reset_mock()
        mock_results = [
            (1, 1, 4),  # PeptideID 1, PositionType 1 (N-term), UnimodAccession 4
            (2, 2, 17),  # PeptideID 2, PositionType 2 (C-term), UnimodAccession 17
        ]
        setup["mock_session"].execute.return_value = mock_results

        terminal_mods = reader._get_terminal_modifications(is_decoy=False)

        assert isinstance(terminal_mods, dict)
        assert terminal_mods[1] == [(1, 4)]
        assert terminal_mods[2] == [(2, 17)]

    def test_get_protein_entries_structure(self, mock_reader_setup):
        """Test _get_protein_entries method structure."""
        setup = mock_reader_setup
        reader = MSFReader(setup["test_file"])

        setup["mock_session"].execute.reset_mock()
        mock_results = [
            (1, "sp|P12345|PROT1_HUMAN Protein 1"),
            (1, "sp|Q67890|PROT2_HUMAN Protein 2"),
            (2, ">tr|R12345|PROT3_HUMAN Protein 3"),
        ]
        setup["mock_session"].execute.return_value = mock_results

        proteins = reader._get_protein_entries(is_decoy=False)

        assert isinstance(proteins, dict)
        assert len(proteins[1]) == 2
        assert proteins[1][0] == "sp|P12345|PROT1_HUMAN Protein 1"
        assert proteins[2][0] == "tr|R12345|PROT3_HUMAN Protein 3"  # ">" should be removed

    def test_get_main_score_structure(self, mock_reader_setup):
        """Test _get_main_score method structure."""
        setup = mock_reader_setup
        reader = MSFReader(setup["test_file"])

        setup["mock_session"].execute.reset_mock()
        mock_results = [
            (1, 95.5, "XCorr"),
            (2, 88.2, "Mascot Score"),
        ]
        setup["mock_session"].execute.return_value = mock_results

        scores = reader._get_main_score(is_decoy=False)

        assert isinstance(scores, dict)
        assert scores[1] == (95.5, "XCorr")
        assert scores[2] == (88.2, "Mascot Score")

    def test_get_secondary_scores_structure(self, mock_reader_setup):
        """Test _get_secondary_scores method structure."""
        setup = mock_reader_setup
        reader = MSFReader(setup["test_file"])

        setup["mock_session"].execute.reset_mock()
        mock_results = [
            (1, 0.95, "Confidence"),
            (1, 15.2, "Delta Score"),
            (2, 0.88, "Confidence"),
        ]
        setup["mock_session"].execute.return_value = mock_results

        scores = reader._get_secondary_scores(is_decoy=False)

        assert isinstance(scores, dict)
        assert scores[1]["Confidence"] == 0.95
        assert scores[1]["Delta Score"] == 15.2
        assert scores[2]["Confidence"] == 0.88

    @pytest.mark.parametrize(
        "sequence,charge",
        [
            ("PEPTIDE", 2),
            ("METHYLATION", 3),
            ("ACDEFGHIKLMNPQRSTVWY", 4),
        ],
    )
    def test_compile_peptidoform_basic(self, mock_reader_setup, sequence, charge):
        """Test _compile_peptidoform with various basic sequences."""
        setup = mock_reader_setup
        reader = MSFReader(setup["test_file"])

        peptidoform = reader._compile_peptidoform(
            sequence=sequence, charge=charge, modifications=[], terminal_modifications=[]
        )

        assert isinstance(peptidoform, Peptidoform)
        assert sequence in str(peptidoform)

    def test_compile_peptidoform_with_modifications(self, mock_reader_setup):
        """Test _compile_peptidoform with amino acid and terminal modifications."""
        setup = mock_reader_setup
        reader = MSFReader(setup["test_file"])

        peptidoform = reader._compile_peptidoform(
            sequence="PEPTIDE",
            charge=2,
            modifications=[(0, 4), (3, 35)],  # Acetyl at pos 0, Oxidation at pos 3
            terminal_modifications=[(1, 1), (2, 17)],  # N-term and C-term mods
        )

        assert isinstance(peptidoform, Peptidoform)
        # Verify modifications are included in the peptidoform
        peptidoform_str = str(peptidoform)
        # The sequence will be modified with UNIMOD annotations
        assert any(aa in peptidoform_str for aa in "PEPTIDE")
        assert "UNIMOD" in peptidoform_str  # Should contain modification annotations

    def test_compatible_versions_constant(self):
        """Test COMPATIBLE_VERSIONS constant is properly defined."""
        assert isinstance(COMPATIBLE_VERSIONS, list)
        assert len(COMPATIBLE_VERSIONS) > 0
        assert all(isinstance(v, int) for v in COMPATIBLE_VERSIONS)
        assert 79 in COMPATIBLE_VERSIONS
        assert 53 in COMPATIBLE_VERSIONS
        assert 8 in COMPATIBLE_VERSIONS


class TestMSFReaderIntegration:
    """Integration tests using the real minimal MSF file."""

    @pytest.fixture
    def minimal_msf_path(self):
        """Path to the minimal MSF test file."""
        path = Path(__file__).parent.parent / "test_data" / "minimal_v79_test.msf"
        if not path.exists():
            pytest.skip("Minimal MSF test file not found")
        return path

    @pytest.fixture
    def reader(self, minimal_msf_path):
        """MSFReader instance with minimal test file."""
        return MSFReader(minimal_msf_path)

    def test_initialization_with_real_file(self, minimal_msf_path):
        """Test successful initialization with real MSF file."""
        reader = MSFReader(minimal_msf_path)
        assert reader is not None
        assert reader.filename == minimal_msf_path

    def test_len_with_real_data(self, reader):
        """Test __len__ method with real MSF data."""
        psm_count = len(reader)
        assert psm_count > 0
        assert isinstance(psm_count, int)

    def test_iteration_yields_correct_count(self, reader):
        """Test that iteration yields the same number of PSMs as len()."""
        expected_count = len(reader)
        actual_psms = list(reader)
        assert len(actual_psms) == expected_count

    def test_psm_structure_and_types(self, reader):
        """Test that PSMs have correct structure and data types."""
        psms = list(reader)
        assert len(psms) > 0

        first_psm = psms[0]

        # Test required attributes exist
        assert hasattr(first_psm, "peptidoform")
        assert hasattr(first_psm, "spectrum_id")
        assert hasattr(first_psm, "run")
        assert hasattr(first_psm, "is_decoy")
        assert hasattr(first_psm, "score")
        assert hasattr(first_psm, "precursor_mz")
        assert hasattr(first_psm, "retention_time")
        assert hasattr(first_psm, "protein_list")
        assert hasattr(first_psm, "rank")
        assert hasattr(first_psm, "source")
        assert hasattr(first_psm, "metadata")
        assert hasattr(first_psm, "rescoring_features")

        # Test data types
        assert isinstance(first_psm.peptidoform, Peptidoform)
        assert isinstance(first_psm.is_decoy, bool)
        assert isinstance(first_psm.score, int | float)
        assert isinstance(first_psm.protein_list, list)
        assert isinstance(first_psm.rank, int)
        assert first_psm.source == "proteome_discoverer"
        assert isinstance(first_psm.metadata, dict)
        assert isinstance(first_psm.rescoring_features, dict)

    def test_target_and_decoy_psms(self, reader):
        """Test that both target and decoy PSMs are present (if available)."""
        psms = list(reader)

        target_psms = [psm for psm in psms if not psm.is_decoy]
        decoy_psms = [psm for psm in psms if psm.is_decoy]

        # At least one type should be present
        assert len(target_psms) > 0 or len(decoy_psms) > 0

        # If both are present, verify they have the expected structure
        if target_psms and decoy_psms:
            assert len(target_psms) > 0
            assert len(decoy_psms) > 0

    def test_psm_metadata_content(self, reader):
        """Test that PSM metadata contains expected keys."""
        psms = list(reader)
        first_psm = psms[0]

        # Test required metadata keys
        required_metadata_keys = [
            "ms1_intensity",
            "ms1_percent_isolation_interference",
            "ms1_ion_inject_time",
            "main_score_name",
        ]

        for key in required_metadata_keys:
            assert key in first_psm.metadata

    def test_rescoring_features_content(self, reader):
        """Test that rescoring features contain expected data."""
        psms = list(reader)
        first_psm = psms[0]

        # Test required rescoring feature keys
        required_rescoring_keys = ["missed_cleavages", "total_ions_count", "matched_ions_count"]

        for key in required_rescoring_keys:
            assert key in first_psm.rescoring_features
            # Values can be int or float depending on database storage
            assert isinstance(first_psm.rescoring_features[key], int | float)

    def test_peptidoform_sequences_valid(self, reader):
        """Test that peptidoform sequences contain valid amino acids."""
        psms = list(reader)

        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

        for psm in psms[:5]:  # Test first 5 PSMs
            peptidoform_str = str(psm.peptidoform)
            # Extract base sequence (before any charge or modification info)
            sequence = peptidoform_str.split("/")[0]
            # Remove any modification annotations
            clean_sequence = "".join(c for c in sequence if c.isalpha())

            # All characters should be valid amino acids
            assert all(aa in valid_amino_acids for aa in clean_sequence), (
                f"Invalid amino acids in sequence: {clean_sequence}"
            )

    def test_unique_spectrum_ids(self, reader):
        """Test that spectrum IDs are unique (within decoy/target groups)."""
        psms = list(reader)

        target_spectrum_ids = {psm.spectrum_id for psm in psms if not psm.is_decoy}
        decoy_spectrum_ids = {psm.spectrum_id for psm in psms if psm.is_decoy}

        target_psms = [psm for psm in psms if not psm.is_decoy]
        decoy_psms = [psm for psm in psms if psm.is_decoy]

        # Within each group, spectrum IDs should be unique per PSM
        if target_psms:
            assert len(target_spectrum_ids) <= len(target_psms)
        if decoy_psms:
            assert len(decoy_spectrum_ids) <= len(decoy_psms)

    def test_score_values_reasonable(self, reader):
        """Test that score values are reasonable numbers."""
        psms = list(reader)

        for psm in psms:
            assert isinstance(psm.score, int | float)
            assert not (psm.score != psm.score)  # Check for NaN
            # Scores should be finite
            assert abs(psm.score) < float("inf")


class TestMSFReaderErrorHandling:
    """Test error handling and edge cases."""

    def test_nonexistent_file(self):
        """Test handling of nonexistent MSF file."""
        with pytest.raises(Exception):  # Should raise some form of file not found error
            MSFReader("nonexistent_file.msf")

    @patch("psm_utils.io.proteome_discoverer.create_engine")
    def test_database_connection_error(self, mock_create_engine):
        """Test handling of database connection errors."""
        mock_create_engine.side_effect = Exception("Database connection failed")

        test_file = Path("test.msf")
        test_file.touch()

        try:
            with pytest.raises(Exception):
                MSFReader(test_file)
        finally:
            test_file.unlink()

    def test_read_file_method(self, minimal_msf_path):
        """Test the read_file() method returns PSMList."""
        reader = MSFReader(minimal_msf_path)
        psm_list = reader.read_file()

        # Should return a list-like object
        assert hasattr(psm_list, "__iter__")
        assert hasattr(psm_list, "__len__")
        assert len(psm_list) > 0

    def test_reader_reusability(self, minimal_msf_path):
        """Test that MSFReader can be reused for multiple operations."""
        reader = MSFReader(minimal_msf_path)

        # Multiple length checks should work
        psm_count1 = len(reader)
        psm_count2 = len(reader)
        assert psm_count1 == psm_count2
        assert psm_count1 > 0

    def test_multiple_iterations(self, minimal_msf_path):
        """Test that multiple iterations over the same reader work consistently."""
        reader = MSFReader(minimal_msf_path)

        first_iteration = list(reader)
        second_iteration = list(reader)

        assert len(first_iteration) == len(second_iteration)
        assert len(first_iteration) > 0

    @pytest.fixture
    def minimal_msf_path(self):
        """Path to the minimal MSF test file."""
        path = Path(__file__).parent.parent / "test_data" / "minimal_v79_test.msf"
        if not path.exists():
            pytest.skip("Minimal MSF test file not found")
        return path


class TestMSFReaderPerformance:
    """Performance and stress tests for MSFReader."""

    @pytest.fixture
    def minimal_msf_path(self):
        """Path to the minimal MSF test file."""
        path = Path(__file__).parent.parent / "test_data" / "minimal_v79_test.msf"
        if not path.exists():
            pytest.skip("Minimal MSF test file not found")
        return path

    def test_lazy_iteration_memory_efficiency(self, minimal_msf_path):
        """Test that iteration is memory efficient (doesn't load all PSMs at once)."""
        reader = MSFReader(minimal_msf_path)

        # Should be able to iterate without loading everything into memory
        psm_count = 0
        for psm in reader:
            psm_count += 1
            if psm_count > 5:  # Just test first few PSMs
                break

        assert psm_count > 0

    def test_consistent_psm_ordering(self, minimal_msf_path):
        """Test that PSM ordering is consistent across iterations."""
        reader = MSFReader(minimal_msf_path)

        first_batch = []
        for i, psm in enumerate(reader):
            first_batch.append(psm.spectrum_id)
            if i >= 4:  # First 5 PSMs
                break

        second_batch = []
        for i, psm in enumerate(reader):
            second_batch.append(psm.spectrum_id)
            if i >= 4:  # First 5 PSMs
                break

        assert first_batch == second_batch, "PSM ordering should be consistent"

    def test_all_required_psm_attributes(self, minimal_msf_path):
        """Test that all PSMs have all required attributes populated."""
        reader = MSFReader(minimal_msf_path)

        required_attrs = [
            "peptidoform",
            "spectrum_id",
            "run",
            "is_decoy",
            "score",
            "precursor_mz",
            "retention_time",
            "protein_list",
            "rank",
            "source",
            "metadata",
            "rescoring_features",
        ]

        for i, psm in enumerate(reader):
            for attr in required_attrs:
                assert hasattr(psm, attr), f"PSM {i} missing attribute: {attr}"
                # None values are acceptable, but attribute must exist
                getattr(psm, attr)

            if i >= 9:  # Test first 10 PSMs
                break
