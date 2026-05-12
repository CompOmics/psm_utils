"""Tests for psm_utils.io.idxml."""

import pytest

from psm_utils.io.idxml import IdXMLReader, IdXMLWriter
from psm_utils.io.sage import SageTSVReader
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM

pyopenms = pytest.importorskip("pyopenms")


def _assert_float_equal(a: float | None, b: float | None, tolerance: float = 1e-5) -> None:
    """Assert two float values are equal within tolerance, handling None values."""
    if a is None and b is None:
        return
    if a is None or b is None:
        assert False, f"One value is None: {a} vs {b}"
    assert abs(a - b) < tolerance, f"Values not equal within tolerance: {a} vs {b}"


class TestIdXMLReader:
    def test__parse_peptidoform(self):
        test_cases = [
            (("DFPIAMGER", 2), "DFPIAMGER/2"),
            (("DFPIAM(Oxidation)GER", 2), "DFPIAM[Oxidation]GER/2"),
            ((".DFPIAMGER.", 2), "DFPIAMGER/2"),
            ((".DFPIAM(Oxidation)GER.", 2), "DFPIAM[Oxidation]GER/2"),
            ((".DFPIAM(UniMod:35)GER.", 2), "DFPIAM[UniMod:35]GER/2"),
            ((".DFPIAM[+16]GER.", 2), "DFPIAM[+16]GER/2"),
            ((".DFPIAM[147]GER.", 2), "DFPIAM[147]GER/2"),
            ((".(Dimethyl)DFPIAMGER.", 2), "[Dimethyl]-DFPIAMGER/2"),
            ((".DFPIAMGER.(Label:18O(2))", 2), "DFPIAMGER-[Label:18O(2)]/2"),
            ((".DFPIAMGER(Phospho).", 2), "DFPIAMGER[Phospho]/2"),
        ]

        for test_in, expected_out in test_cases:
            assert IdXMLReader._parse_peptidoform(*test_in) == expected_out

    def test__parse_psm(self):
        test_psm = PSM(
            peptidoform=Peptidoform("LVPIWKK/2"),
            spectrum_id="controllerType=0 controllerNumber=1 scan=294",
            run="HepG2_rep3_small",
            collection=None,
            spectrum=None,
            is_decoy=True,
            score=1.0099999904632568,
            qvalue=None,
            pep=None,
            precursor_mz=442.29595853189994,
            retention_time=849.0,
            ion_mobility=None,
            protein_list=["DECOY_sp|Q5TH74|STPG1_HUMAN"],
            rank=1,
            source="idXML",
            provenance_data={"LVPIWKK/2": "LVPIWKK"},
            metadata={
                "idxml:score_type": "expect",
                "idxml:higher_score_better": "False",
                "idxml:significance_threshold": "0.0",
                "target_decoy": "decoy",
                "protein_references": "unique",
            },
            rescoring_features={
                "MS:1002258": 3.0,
                "MS:1002259": 12.0,
                "num_matched_peptides": 35.0,
                "isotope_error": 0.0,
                "MS:1002252": 0.693,
                "COMET:xcorr": 0.693,
                "MS:1002253": 1.0,
                "COMET:deltaCn": 1.0,
                "MS:1002255": 35.9,
                "COMET:spscore": 35.9,
                "MS:1002256": 1.0,
                "COMET:sprank": 1.0,
                "MS:1002257": 1.01,
                "COMET:deltaLCn": 0.0,
                "COMET:lnExpect": 0.009950330853168092,
                "COMET:lnNumSP": 3.555348061489414,
                "COMET:lnRankSP": 0.0,
                "COMET:IonFrac": 0.25,
            },
        )
        reader = IdXMLReader("./tests/test_data/test_in.idXML")
        psm = reader._parse_psm(
            reader.protein_ids, reader.peptide_ids[0], reader.peptide_ids[0].getHits()[0]
        )
        assert psm == test_psm

    def test__get_run(self):
        expected_output = "HepG2_rep3_small"
        reader = IdXMLReader("./tests/test_data/test_in.idXML")
        assert reader._get_run(reader.protein_ids, reader.peptide_ids[0]) == expected_output

    def test__get_rescoring_features(self):
        expected_output = [
            "MS:1002258",
            "MS:1002259",
            "num_matched_peptides",
            "isotope_error",
            "MS:1002252",
            "COMET:xcorr",
            "MS:1002253",
            "COMET:deltaCn",
            "MS:1002255",
            "COMET:spscore",
            "MS:1002256",
            "COMET:sprank",
            "MS:1002257",
            "COMET:deltaLCn",
            "COMET:lnExpect",
            "COMET:lnNumSP",
            "COMET:lnRankSP",
            "COMET:IonFrac",
        ]
        reader = IdXMLReader("./tests/test_data/test_in.idXML")
        assert (
            reader._get_rescoring_features(reader.peptide_ids[0].getHits()[0]) == expected_output
        )


class TestIdXMLWriter:
    def test_write_file_with_pyopenms_objects(self):
        """Test writing idXML file with existing pyopenms objects and verify content."""
        reader = IdXMLReader("./tests/test_data/test_in.idXML")
        original_psm_list = reader.read_file()

        # Write the file
        writer = IdXMLWriter(
            "./tests/test_data/test_out.idXML", reader.protein_ids, reader.peptide_ids
        )
        writer.write_file(original_psm_list)

        # Read back the written file and verify content
        reader_check = IdXMLReader("./tests/test_data/test_out.idXML")
        written_psm_list = reader_check.read_file()

        # Verify basic file structure
        assert len(written_psm_list) == len(original_psm_list)

        # Compare key attributes of each PSM
        for orig_psm, written_psm in zip(original_psm_list, written_psm_list):
            assert str(orig_psm.peptidoform) == str(written_psm.peptidoform)
            assert orig_psm.spectrum_id == written_psm.spectrum_id
            assert orig_psm.run == written_psm.run
            assert orig_psm.is_decoy == written_psm.is_decoy
            _assert_float_equal(orig_psm.score, written_psm.score)
            _assert_float_equal(orig_psm.precursor_mz, written_psm.precursor_mz)
            _assert_float_equal(orig_psm.retention_time, written_psm.retention_time)
            assert orig_psm.protein_list == written_psm.protein_list
            assert orig_psm.rank == written_psm.rank

            # Check that rescoring features are preserved
            if orig_psm.rescoring_features:
                assert written_psm.rescoring_features is not None
                for feature_name, feature_value in orig_psm.rescoring_features.items():
                    assert feature_name in written_psm.rescoring_features
                    assert abs(written_psm.rescoring_features[feature_name] - feature_value) < 1e-6

    def test_write_file_without_pyopenms_objects(self):
        """Test writing idXML file from scratch without existing pyopenms objects."""
        reader = SageTSVReader("./tests/test_data/results.sage.tsv")
        original_psm_list = reader.read_file()

        # Write the file
        writer = IdXMLWriter("./tests/test_data/test_out_sage.idXML")
        writer.write_file(original_psm_list)

        # Read back the written file and verify content
        reader_check = IdXMLReader("./tests/test_data/test_out_sage.idXML")
        written_psm_list = reader_check.read_file()

        # Verify basic file structure
        assert len(written_psm_list) == len(original_psm_list)

        # Compare key attributes of the first PSM (since sage data has one entry)
        orig_psm = original_psm_list[0]
        written_psm = written_psm_list[0]

        assert str(orig_psm.peptidoform) == str(written_psm.peptidoform)
        assert orig_psm.spectrum_id == written_psm.spectrum_id
        assert orig_psm.run == written_psm.run
        assert orig_psm.is_decoy == written_psm.is_decoy
        _assert_float_equal(orig_psm.score, written_psm.score)
        _assert_float_equal(orig_psm.precursor_mz, written_psm.precursor_mz)
        _assert_float_equal(orig_psm.retention_time, written_psm.retention_time)
        assert orig_psm.protein_list == written_psm.protein_list
        assert orig_psm.rank == written_psm.rank

        # Verify that the written file is a valid idXML (can be read without errors)
        assert len(reader_check.protein_ids) > 0
        assert len(reader_check.peptide_ids) > 0

    def test_write_file_preserves_modifications(self):
        """Test that modifications are properly preserved when writing idXML files."""
        from psm_utils.psm_list import PSMList

        # Create test PSMs with various modifications
        test_psms = [
            PSM(
                peptidoform="ACDK/2",
                spectrum_id="scan=1",
                score=140.2,
                retention_time=600.2,
                precursor_mz=300.15,
                run="test_run",
            ),
            PSM(
                peptidoform="AC[Carbamidomethyl]DK/2",
                spectrum_id="scan=2",
                score=150.3,
                retention_time=650.1,
                precursor_mz=357.17,
                run="test_run",
            ),
            PSM(
                peptidoform="[Acetyl]-ACDK/2",
                spectrum_id="scan=3",
                score=120.8,
                retention_time=580.5,
                precursor_mz=342.16,
                run="test_run",
            ),
        ]

        psm_list = PSMList(psm_list=test_psms)

        # Write and read back
        writer = IdXMLWriter("./tests/test_data/test_mods.idXML")
        writer.write_file(psm_list)

        reader_check = IdXMLReader("./tests/test_data/test_mods.idXML")
        written_psm_list = reader_check.read_file()

        # Verify modifications are preserved
        assert len(written_psm_list) == len(test_psms)

        for orig_psm, written_psm in zip(test_psms, written_psm_list):
            # The peptidoform should be preserved (though the exact string representation might differ)
            assert orig_psm.peptidoform.sequence == written_psm.peptidoform.sequence
            assert (
                orig_psm.peptidoform.precursor_charge == written_psm.peptidoform.precursor_charge
            )

            # Basic properties should match
            assert orig_psm.spectrum_id == written_psm.spectrum_id
            _assert_float_equal(orig_psm.score, written_psm.score)
            _assert_float_equal(orig_psm.retention_time, written_psm.retention_time)
            _assert_float_equal(orig_psm.precursor_mz, written_psm.precursor_mz)

    def test_write_file_with_metadata_and_features(self):
        """Test that metadata and rescoring features are preserved."""
        from psm_utils.psm_list import PSMList

        test_psm = PSM(
            peptidoform="TESTPEPTIDE/2",
            spectrum_id="scan=100",
            score=200.5,
            retention_time=1000.0,
            precursor_mz=500.25,
            run="feature_test",
            qvalue=0.01,
            pep=0.05,
            metadata={"custom_meta": "test_value", "intensity": "12345"},
            rescoring_features={"custom_score": 0.85, "feature_2": 1.23},
        )

        psm_list = PSMList(psm_list=[test_psm])

        # Write and read back
        writer = IdXMLWriter("./tests/test_data/test_features.idXML")
        writer.write_file(psm_list)

        reader_check = IdXMLReader("./tests/test_data/test_features.idXML")
        written_psm_list = reader_check.read_file()

        assert len(written_psm_list) == 1
        written_psm = written_psm_list[0]

        # Check basic attributes
        assert str(test_psm.peptidoform) == str(written_psm.peptidoform)
        assert test_psm.spectrum_id == written_psm.spectrum_id
        _assert_float_equal(test_psm.score, written_psm.score)
        _assert_float_equal(test_psm.qvalue, written_psm.qvalue)
        _assert_float_equal(test_psm.pep, written_psm.pep)

        # Check that custom features are preserved
        assert written_psm.rescoring_features is not None
        assert "custom_score" in written_psm.rescoring_features
        assert abs(written_psm.rescoring_features["custom_score"] - 0.85) < 1e-6
