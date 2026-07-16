"""Tests for psm_utils.io.diann."""

from psm_utils.io.diann import DIANNParquetReader, DIANNTSVReader
from psm_utils.psm import PSM

test_psm = PSM(
    peptidoform="AAAAEINVKDPKEDLETSVVDEGR/4",
    spectrum_id="116903",
    run="LFQ_Orbitrap_AIF_Yeast_03",
    collection=None,
    spectrum=None,
    is_decoy=False,
    score=0.995107,
    qvalue=0.000548193,
    pep=0.0104343,
    precursor_mz=None,
    retention_time=75.2574,
    ion_mobility=0,
    protein_list=["P38156"],
    rank=None,
    source="diann",
    metadata={},
    rescoring_features={
        "RT": 75.2574,
        "Predicted.RT": 75.2713,
        "iRT": 33.9222,
        "Predicted.iRT": 33.8999,
        "Ms1.Profile.Corr": 0.347567,
        "Ms1.Area": 849940,
        "IM": 0,
        "iIM": 0,
        "Predicted.IM": 0,
        "Predicted.iIM": 0,
    },
)


class TestDIANNTSVReader:
    def test_iter(self):
        with DIANNTSVReader("./tests/test_data/test_diann.tsv") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_psm

    def test__parse_peptidoform(self):
        test_cases = [
            (("ACDE", "4"), "ACDE/4"),
            (("AC(UniMod:1)DE", "4"), "AC[UNIMOD:1]DE/4"),
            (("(UniMod:4)ACDE", "4"), "[UNIMOD:4]-ACDE/4"),
        ]

        reader = DIANNTSVReader("./tests/test_data/test_diann.tsv")
        for (peptide, charge), expected in test_cases:
            assert reader._parse_peptidoform(peptide, charge) == expected


test_parquet_psm = PSM(
    peptidoform="[UNIMOD:1]-AAAAEINVKDPKEDLETSVVDEGR/4",
    spectrum_id="(UniMod:1)AAAAEINVKDPKEDLETSVVDEGR4",
    run="LFQ_Orbitrap_AIF_Yeast_03",
    collection=None,
    spectrum=None,
    is_decoy=False,
    score=None,
    qvalue=0.000548193,
    pep=0.0104343,
    precursor_mz=621.5723,
    retention_time=75.2574,
    ion_mobility=0.0,
    protein_list=["P38156"],
    rank=None,
    source="diann",
    metadata={},
    rescoring_features={
        "RT": 75.2574,
        "Predicted.RT": 75.2713,
        "iRT": 33.9222,
        "Predicted.iRT": 33.8999,
        "Ms1.Profile.Corr": 0.347567,
        "Ms1.Area": 849940.0,
        "IM": 0.0,
        "iIM": 0.0,
        "Predicted.IM": 0.0,
        "Predicted.iIM": 0.0,
    },
)


class TestDIANNParquetReader:
    def test_iter(self):
        with DIANNParquetReader("./tests/test_data/test_diann.parquet") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_parquet_psm
