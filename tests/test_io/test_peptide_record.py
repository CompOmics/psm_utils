"""Tests for psm_utils.io.peptide_record."""

import pytest

from psm_utils.io.peptide_record import (
    InvalidPeprecError,
    InvalidPeprecModificationError,
    _analyze_peprec_file,
    peprec_to_proforma,
)


class TestPeprecFileAnalysis:
    def test_analyze_peprec_file(self):
        # Tab
        separator, header = _analyze_peprec_file("./tests/test_data/peprec.tsv")
        assert separator == "\t"
        assert "spec_id" in header
        assert "peptide" in header
        assert "modifications" in header

        # Comma
        separator, header = _analyze_peprec_file("./tests/test_data/peprec.csv")
        assert separator == ","
        assert "spec_id" in header
        assert "peptide" in header
        assert "modifications" in header

        # Space
        separator, header = _analyze_peprec_file("./tests/test_data/peprec.txt")
        assert separator == " "
        assert "spec_id" in header
        assert "peptide" in header
        assert "modifications" in header

        # Invalid: Mixed use of separators
        with pytest.raises(InvalidPeprecError):
            _analyze_peprec_file("./tests/test_data/peprec_invalid.csv")


class TestPeprecToProforma:
    def test_peprec_to_proforma(self):
        # Valid cases
        valid_test_cases = [
            (("ACDMEK", ""), "ACDMEK"),
            (("ACDMEK", "-"), "ACDMEK"),
            (("ACDMEK", "4|Oxidation"), "ACDM[Oxidation]EK"),
            (("ACDMEK", "0|Acetylation"), "[Acetylation]-ACDMEK"),
            (("ACDMEK", "-1|Amidation"), "ACDMEK-[Amidation]"),
            (("ACDMEK", "4|Oxidation|-1|Amidation"), "ACDM[Oxidation]EK-[Amidation]"),
            (("ACDMEK", "0|Acetyl|4|Ox|-1|Amide"), "[Acetyl]-ACDM[Ox]EK-[Amide]"),
            (
                ("ACDMEK", "6|Methylation|-1|Amide"),
                "ACDMEK[Methylation]-[Amide]",
            ),  # See levitsky/pyteomics/#77
            (("MCDMEK", "0|Acetyl|1|Oxidation"), "[Acetyl]-M[Oxidation]CDMEK"),
            (("ACDMEK", "", "2"), "ACDMEK/2"),
        ]
        for test_in, expected_out in valid_test_cases:
            test_out = peprec_to_proforma(*test_in)
            assert test_out.proforma == expected_out

        # Invalid cases
        invalid_test_cases = [
            ("ACDE", "|"),
            ("ACDE", "8|Oxidation"),
        ]
        for test_in in invalid_test_cases:
            with pytest.raises(InvalidPeprecModificationError):
                peprec_to_proforma(*test_in)


class TestPeptideRecordReader:
    # TODO!!!
    def test_read_file(self):
        pass
