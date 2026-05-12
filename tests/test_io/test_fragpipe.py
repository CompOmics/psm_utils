"""Tests for psm_utils.io.fragpipe."""

from psm_utils.io.fragpipe import FragPipeReader
from psm_utils.psm import PSM

test_psm = PSM(
    peptidoform="TGAPNNGQYGADNGNPNGER/2",
    spectrum_id="00001",
    run="LFQ_Orbitrap_AIF_Yeast_01_Q1",
    collection=None,
    spectrum=None,
    is_decoy=False,
    score=57.2940,
    qvalue=None,
    pep=1 - 1.0000,
    precursor_mz=1001.9336,
    retention_time=2432.1640,
    ion_mobility=None,
    protein_list=["sp|P40159|YNU8_YEAST"],
    rank=None,
    source="FragPipe",
    metadata={},
    rescoring_features={
        "Peptide Length": 20,
        "Retention": 2432.1640,
        "Observed Mass": 2001.8539,
        "Observed M/Z": 1001.9342,
        "Calculated Peptide Mass": 2001.8524,
        "Calculated M/Z": 1001.9335,
        "Delta Mass": 0.0002,
        "Number of Missed Cleavages": 0,
    },
)


class TestFragpipeReader:
    def test_iter(self):
        with FragPipeReader("./tests/test_data/test_fragpipe.tsv") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_psm

    def test__parse_peptidoform(self):
        test_cases = [
            (("LHMTNQNMEK", "3M(15.994915), C-term(17.034480)", "3"), "LHM[+15.9949]TNQNMEK-[+17.0345]/3"),
            (("ANIAVQR", "N-term(42.0106)", "2"), "[+42.0106]-ANIAVQR/2"),
            (("IPAVTYPK", "", "2"), "IPAVTYPK/2"),
        ]

        reader = FragPipeReader("./tests/test_data/test_fragpipe.tsv")
        for (peptide, modifications, charge), expected in test_cases:
            assert reader._parse_peptidoform(peptide, modifications, charge) == expected

    def test__parse_spectrum_id(self):
        test_cases = [
            ("LFQ_Orbitrap_AIF_Human_01.101124.101124.0", "101124"),
            ("LFQ.Orbitrap.AIF.Human.01.101124.101124.0", "101124"),
            ("101124", "101124"),
        ]

        reader = FragPipeReader("./tests/test_data/test_fragpipe.tsv")
        for spectrum, expected in test_cases:
            assert reader._parse_spectrum_id(spectrum) == expected
