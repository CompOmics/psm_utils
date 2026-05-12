"""Tests for psm_utils.io.percolator."""

from unittest.mock import mock_open, patch

import pytest

from psm_utils.io.percolator import PercolatorIOException, PercolatorTabReader, PercolatorTabWriter


class TestPercolatorTabReader:
    def test__infer_charge_columns(self):
        test_cases = [
            (["Charge"], ("Charge", {})),
            (["charge"], ("charge", {})),
            (["Charge", "Charge2"], ("Charge", {2: "Charge2"})),
            (
                ["Charge2", "Charge3", "Charge4"],
                (None, {2: "Charge2", 3: "Charge3", 4: "Charge4"}),
            ),
            (
                ["charge2", "charge3", "charge4"],
                (None, {2: "charge2", 3: "charge3", 4: "charge4"}),
            ),
        ]

        for test_in, expected_out in test_cases:
            assert expected_out == PercolatorTabReader._infer_charge_columns(test_in)

    def test_parse_peptidoform(self):
        test_cases = [
            # Basic cases
            (("ACDEFGHR", None), "ACDEFGHR"),
            (("K.ACDEFGHR.I", 1), "ACDEFGHR/1"),
            (("K.ACDEFGHR.-", 2), "ACDEFGHR/2"),
            (("-.ACDEFGHR.I", 3), "ACDEFGHR/3"),
            (("-.ACDEFGHR.-", None), "ACDEFGHR"),
            # N-terminal modifications
            (("-.n[42.0106]ACDEFGHR.-", None), "[+42.0106]-ACDEFGHR"),
            (("n[42.0106]ACDEFGHR", None), "[+42.0106]-ACDEFGHR"),  # Without flanking
            (("-.n[43]ACDEFGHR.-", 2), "[+43]-ACDEFGHR/2"),  # Integer mass
            # C-terminal modifications
            (("-.ACDEFGHRc[-0.9840].-", None), "ACDEFGHR-[-0.984]"),
            (("ACDEFGHRc[-0.9840]", None), "ACDEFGHR-[-0.984]"),  # Without flanking
            (("-.ACDEFGHRc[17.0265].-", 2), "ACDEFGHR-[+17.0265]/2"),  # Positive C-term
            # Internal modifications
            (("-.ACDEFM[15.9949]GHR.-", None), "ACDEFM[+15.9949]GHR"),
            (("-.ACDEM[-18.010565]GHR.-", None), "ACDEM[-18.010565]GHR"),  # Negative internal
            (("-.AC[57.021]DEFGHR.-", None), "AC[+57.021]DEFGHR"),  # Carbamidomethyl
            # Multiple modifications
            (("-.n[43]ACDEFM[16]GHR.-", None), "[+43]-ACDEFM[+16]GHR"),  # N-term + internal
            (("-.ACDEFM[16]GHRc[-1].-", None), "ACDEFM[+16]GHR-[-1]"),  # Internal + C-term
            (("-.n[42]ACDEFM[16]GHRc[-1].-", 2), "[+42]-ACDEFM[+16]GHR-[-1]/2"),  # All three
            (("-.AC[57]DEM[16]GHK.-", None), "AC[+57]DEM[+16]GHK"),  # Multiple internal
            # Already has '+' sign (should not add another)
            (("-.ACDEFM[+15.9949]GHR.-", None), "ACDEFM[+15.9949]GHR"),
            (("-.n[+42.0106]ACDEFGHR.-", None), "[+42.0106]-ACDEFGHR"),
        ]
        for test_in, expected_out in test_cases:
            assert expected_out == PercolatorTabReader._parse_peptidoform(*test_in).proforma


@pytest.fixture
def valid_pout_file():
    return (
        "PSMId\tLabel\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n"
        "PSM1\t1\t0.8\t0.01\t0.001\tPEPTIDE\tPROT1\n"
        "PSM2\t-1\t0.7\t0.02\t0.002\tPEPTIDER\tPROT2\n"
    )


@pytest.fixture
def invalid_pout_file():
    return "Not a valid header\nSome invalid content"


@pytest.fixture
def empty_file():
    return ""


def test_parse_existing_file_valid(valid_pout_file):
    with patch("builtins.open", mock_open(read_data=valid_pout_file)):
        fieldnames, last_scannr = PercolatorTabWriter._parse_existing_file("dummy_path", "pout")
        assert fieldnames == [
            "PSMId",
            "Label",
            "score",
            "q-value",
            "posterior_error_prob",
            "peptide",
            "proteinIds",
        ]
        assert last_scannr == -1  # No ScanNr in POUT style, so it should be -1 by default


def test_parse_existing_file_invalid(invalid_pout_file):
    with patch("builtins.open", mock_open(read_data=invalid_pout_file)):
        with pytest.raises(PercolatorIOException):
            PercolatorTabWriter._parse_existing_file("dummy_path", "pin")


def test_parse_existing_file_empty(empty_file):
    with patch("builtins.open", mock_open(read_data=empty_file)):
        with pytest.raises(PercolatorIOException):
            PercolatorTabWriter._parse_existing_file("dummy_path", "pin")


def test_parse_existing_file_no_scannr_column():
    # Simulate a file without ScanNr column but a valid header and entries
    data = (
        "PSMId\tLabel\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n"
        "PSM1\t1\t0.8\t0.01\t0.001\tPEPTIDE\tPROT1\n"
    )
    with patch("builtins.open", mock_open(read_data=data)):
        fieldnames, last_scannr = PercolatorTabWriter._parse_existing_file("dummy_path", "pout")
        assert fieldnames == [
            "PSMId",
            "Label",
            "score",
            "q-value",
            "posterior_error_prob",
            "peptide",
            "proteinIds",
        ]
        assert last_scannr == -1  # Since no ScanNr column exists, should default to -1


def test_parse_existing_file_with_scannr():
    # Simulate a file with ScanNr column and entries
    data = (
        "PSMId\tLabel\tScanNr\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n"
        "PSM1\t1\t0\t0.8\t0.01\t0.001\tPEPTIDE\tPROT1\n"
        "PSM2\t1\t1\t0.7\t0.02\t0.002\tPEPTIDER\tPROT2\n"
    )
    with patch("builtins.open", mock_open(read_data=data)):
        fieldnames, last_scannr = PercolatorTabWriter._parse_existing_file("dummy_path", "pout")
        assert fieldnames == [
            "PSMId",
            "Label",
            "ScanNr",
            "score",
            "q-value",
            "posterior_error_prob",
            "peptide",
            "proteinIds",
        ]
        assert last_scannr == 1  # The last ScanNr should be 1
