"""Tests for psm_utils.io.json."""

import tempfile
from pathlib import Path

from psm_utils.io.json import JSONReader, JSONWriter
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class TestJSONReader:
    def test__parse_entry(self):
        """Test parsing a JSON entry."""
        entry = {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
        }
        expected = {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
            "provenance_data": {},
            "metadata": {},
            "rescoring_features": {},
        }
        assert JSONReader._parse_entry(entry) == expected

    def test__parse_entry_with_dicts(self):
        """Test parsing a JSON entry with dict attributes."""
        entry = {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
            "provenance_data": {"test": "value"},
            "metadata": {"key": "value"},
            "rescoring_features": {"score": "1.5"},
        }
        expected = entry.copy()
        assert JSONReader._parse_entry(entry) == expected

    def test_iter(self):
        """Test iterating over JSON reader."""
        # Create a temporary JSON file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            f.write('[{"peptidoform": "ACDEK/2", "spectrum_id": "peptide1"}]')
            temp_path = f.name

        try:
            reader = JSONReader(temp_path)
            psms = list(reader)
            assert len(psms) == 1
            assert psms[0].peptidoform == Peptidoform("ACDEK/2")
            assert psms[0].spectrum_id == "peptide1"
        finally:
            Path(temp_path).unlink()


class TestJSONWriter:
    def test_write_psm(self):
        """Test writing a PSM to JSON."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = f.name

        try:
            psm = PSM(
                peptidoform="ACDEK/2",
                spectrum_id="peptide1",
                run="test_run",
                score=1.5,
            )

            with JSONWriter(temp_path, indent=2) as writer:
                writer.write_psm(psm)

            # Read back and verify
            with JSONReader(temp_path) as reader:
                psms = list(reader)
                assert len(psms) == 1
                assert psms[0].peptidoform == Peptidoform("ACDEK/2")
                assert psms[0].spectrum_id == "peptide1"
                assert psms[0].run == "test_run"
                assert psms[0].score == 1.5
        finally:
            Path(temp_path).unlink()

    def test_write_file(self):
        """Test writing a PSMList to JSON."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = f.name

        try:
            psm_list = PSMList(
                psm_list=[
                    PSM(peptidoform="ACDEK/2", spectrum_id="peptide1"),
                    PSM(peptidoform="DEFGH/3", spectrum_id="peptide2"),
                ]
            )

            writer = JSONWriter(temp_path, indent=None)
            writer.write_file(psm_list)

            # Read back and verify
            with JSONReader(temp_path) as reader:
                psms = list(reader)
                assert len(psms) == 2
                assert psms[0].peptidoform == Peptidoform("ACDEK/2")
                assert psms[1].peptidoform == Peptidoform("DEFGH/3")
        finally:
            Path(temp_path).unlink()

    def test_roundtrip(self):
        """Test writing and reading back complex PSMs."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = f.name

        try:
            original_psm = PSM(
                peptidoform="AC[Carbamidomethyl]DEK/2",
                spectrum_id="scan=123",
                run="test_run",
                collection="PXD000001",
                is_decoy=False,
                score=10.5,
                qvalue=0.01,
                pep=0.001,
                precursor_mz=500.5,
                retention_time=100.0,
                ion_mobility=1.5,
                protein_list=["PROT1", "PROT2"],
                rank=1,
                source="test",
                provenance_data={"file": "test.raw"},
                metadata={"instrument": "Orbitrap"},
                rescoring_features={"feature1": 1.0, "feature2": 2.0},
            )

            psm_list = PSMList(psm_list=[original_psm])

            writer = JSONWriter(temp_path, indent=2)
            writer.write_file(psm_list)

            # Read back and verify
            with JSONReader(temp_path) as reader:
                psms = list(reader)
                assert len(psms) == 1
                psm = psms[0]

                assert str(psm.peptidoform) == str(original_psm.peptidoform)
                assert psm.spectrum_id == original_psm.spectrum_id
                assert psm.run == original_psm.run
                assert psm.collection == original_psm.collection
                assert psm.is_decoy == original_psm.is_decoy
                assert psm.score == original_psm.score
                assert psm.qvalue == original_psm.qvalue
                assert psm.pep == original_psm.pep
                assert psm.precursor_mz == original_psm.precursor_mz
                assert psm.retention_time == original_psm.retention_time
                assert psm.ion_mobility == original_psm.ion_mobility
                assert psm.protein_list == original_psm.protein_list
                assert psm.rank == original_psm.rank
                assert psm.source == original_psm.source
                assert psm.provenance_data == original_psm.provenance_data
                assert psm.metadata == original_psm.metadata
        finally:
            Path(temp_path).unlink()
