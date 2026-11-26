import numpy as np
import pytest

from psm_utils import PSM, Peptidoform, PSMList

sample_psm_list = [
    PSM(peptidoform="ACDK", spectrum_id=1, score=140.2),
    PSM(peptidoform="CDEFR", spectrum_id=2, score=132.9),
    PSM(peptidoform="DEM[Oxidation]K", spectrum_id=3, score=55.7),
]


class TestPSMList:
    def test___init__(self):
        _ = PSMList(psm_list=sample_psm_list)

    def test___iter__(self):
        for psm in PSMList(psm_list=sample_psm_list):
            print(psm)

    def test___len__(self):
        psm_list = PSMList(psm_list=sample_psm_list)
        assert len(psm_list) == 3

    def test___get_item__(self):
        psm_list = PSMList(psm_list=sample_psm_list)

        # Single index
        assert psm_list[0] == PSM(peptidoform="ACDK", spectrum_id="1", score=140.2)

        # Slice
        assert psm_list[0:2] == PSMList(
            psm_list=[
                PSM(peptidoform="ACDK", spectrum_id="1", score=140.2),
                PSM(peptidoform="CDEFR", spectrum_id="2", score=132.9),
            ]
        )

        # PSM property as array
        np.testing.assert_equal(psm_list["spectrum_id"], np.array(["1", "2", "3"]))
        np.testing.assert_equal(psm_list["score"], np.array([140.2, 132.9, 55.7]))
        np.testing.assert_equal(psm_list["rank"], np.array([None, None, None]))
        np.testing.assert_equal(psm_list["qvalue"], np.array([np.nan, np.nan, np.nan]))

        # Multiple PSM properties as 2D array
        np.testing.assert_equal(
            psm_list[["spectrum_id", "score"]],
            np.array([["1", 140.2], ["2", 132.9], ["3", 55.7]], dtype=object),
        )

        # Index by multiple indices
        psm_list[0, 2] == PSMList(
            psm_list=[
                PSM(peptidoform="ACDK", spectrum_id="1", score=140.2),
                PSM(peptidoform="DEM[Oxidation]K", spectrum_id="3", score=55.7),
            ]
        )

        # Unsupported indexing type
        with pytest.raises(TypeError):
            psm_list[0.21]

    def test___set_item__(self):
        psm_list = PSMList(psm_list=sample_psm_list)

        # Correct usage
        psm_list["score"] = [10, 20, 30]
        np.testing.assert_equal(psm_list["score"], np.array([10, 20, 30]))

        # Wrong length
        with pytest.raises(ValueError):
            psm_list["score"] = [10, 20, 30, 40]

    def test_set_ranks(self):
        psm_list = PSMList(
            psm_list=[
                PSM(
                    peptidoform=Peptidoform("ACDE/2"),
                    spectrum_id=1,
                    run=None,
                    collection=None,
                    score=3.3,
                ),
                PSM(
                    peptidoform=Peptidoform("ACDF/2"),
                    spectrum_id=2,
                    run=None,
                    collection=None,
                    score=12.0,
                ),
                PSM(
                    peptidoform=Peptidoform("ACDG/2"),
                    spectrum_id=3,
                    run=None,
                    collection=None,
                    score=12.2,
                ),
                PSM(
                    peptidoform=Peptidoform("ACDH/2"),
                    spectrum_id=4,
                    run=None,
                    collection=None,
                    score=2.5,
                ),
            ]
        )

        # Simple unique set of PSMs
        psm_list.set_ranks()
        np.testing.assert_equal(psm_list["rank"], np.array([1, 1, 1, 1]))

        # First two PSMs belong to same spectrum
        psm_list["spectrum_id"] = [1, 1, 2, 3]
        psm_list.set_ranks()
        np.testing.assert_equal(psm_list["rank"], np.array([2, 1, 1, 1]))

        # Single specified run, making spectra unique again
        psm_list["spectrum_id"] = [1, 1, 2, 3]
        psm_list["run"] = [None, "run_1", None, None]
        psm_list.set_ranks()
        np.testing.assert_equal(psm_list["rank"], np.array([1, 1, 1, 1]))

        # All runs specified
        psm_list["spectrum_id"] = [1, 1, 2, 3]
        psm_list["run"] = ["run_1", "run_1", "run_2", "run_2"]
        psm_list.set_ranks()
        np.testing.assert_equal(psm_list["rank"], np.array([2, 1, 1, 1]))

        # Collection makes spectra unique again
        psm_list["spectrum_id"] = [1, 2, 2, 2]
        psm_list["run"] = ["run_1", "run_1", "run_2", "run_2"]
        psm_list["collection"] = ["coll_1", "coll_1", "coll_1", "coll_2"]
        psm_list.set_ranks()
        np.testing.assert_equal(psm_list["rank"], np.array([1, 1, 1, 1]))

        # Lower score is better
        psm_list["spectrum_id"] = [1, 1, 2, 2]
        psm_list["run"] = [None, None, None, None]
        psm_list["collection"] = [None, None, None, None]
        psm_list.set_ranks(lower_score_better=True)
        np.testing.assert_equal(psm_list["rank"], np.array([1, 2, 2, 1]))


def test_rename_modifications_logs_once_for_psm_list(caplog):
    psms = PSMList(
        psm_list=[PSM(peptidoform="M[ox]C[already-custom]K", spectrum_id=str(i)) for i in range(3)]
    )
    with caplog.at_level("INFO", logger="psm_utils.psm_list"):
        psms.rename_modifications({"ox": "Oxidation"})
        psms.rename_modifications({})
    assert [p.peptidoform.proforma for p in psms] == ["M[Oxidation]C[already-custom]K"] * 3
    assert sum("Renaming modifications" in record.message for record in caplog.records) == 1
