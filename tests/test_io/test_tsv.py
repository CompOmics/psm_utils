"""Tests for psm_utils.io.tsv."""

import math

import pytest

from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.io.tsv import TSVReader
from psm_utils.peptidoform import Peptidoform

test_cases = [
    (
        {"peptidoform": "ACDE", "spectrum_id": "1"},
        {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
            "provenance_data": {},
            "metadata": {},
            "rescoring_features": {},
        },
    ),
    (
        {"peptidoform": "ACDE", "spectrum_id": "1", "provenance:test": "value"},
        {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
            "provenance_data": {"test": "value"},
            "metadata": {},
            "rescoring_features": {},
        },
    ),
    (
        # Empty string provenance value (e.g. missing optional provenance field)
        {"peptidoform": "ACDE", "spectrum_id": "1", "provenance:missing": ""},
        {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
            "provenance_data": {"missing": ""},
            "metadata": {},
            "rescoring_features": {},
        },
    ),
    (
        # Empty string rescoring feature value should become NaN
        {"peptidoform": "ACDE", "spectrum_id": "1", "rescoring:score": ""},
        {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
            "provenance_data": {},
            "metadata": {},
            "rescoring_features": {"score": float("nan")},
        },
    ),
]


class TestTSVReader:
    def test__parse_entry(self):
        for test_in, expected_out in test_cases:
            result = TSVReader._parse_entry(test_in)
            for key, expected_val in expected_out.items():
                if isinstance(expected_val, dict):
                    for k, v in expected_val.items():
                        if isinstance(v, float) and math.isnan(v):
                            assert math.isnan(result[key][k])
                        else:
                            assert result[key][k] == v
                else:
                    assert result[key] == expected_val

    def test_iter(self):
        reader = TSVReader("tests/test_data/test.tsv")
        for psm in reader:
            assert psm.peptidoform == Peptidoform("ACDEK/2")
            assert psm.spectrum_id == "peptide1"
            assert psm.provenance_data == {}
            assert psm.metadata == {}
            assert psm.rescoring_features == {}
            break

    def test_iter_raises(self):
        with TSVReader("tests/test_data/peprec.tsv") as reader:
            with pytest.raises(PSMUtilsIOException):
                for psm in reader:
                    pass
