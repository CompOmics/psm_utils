import pytest

import psm_utils.io.maxquant as maxquant

TEST_COL = [
    "Raw file",
    "Scan number",
    "Charge",
    "Modified sequence",
    "Proteins",
    "m/z",
    "Reverse",
    "Retention time",
    "PEP",
    "Score",
]


class TestMSMSReader:
    def test_evaluate_columns(self):
        columns = TEST_COL.copy()
        # Test with the right column names
        maxquant.MSMSReader._evaluate_columns(columns)

        # Test with right columns names but lowercase columnname
        columns[0] = "raw file"
        maxquant.MSMSReader._evaluate_columns(columns)

        # Test when column name is missing
        columns.remove("m/z")
        with pytest.raises(maxquant.MSMSParsingError):
            maxquant.MSMSReader._evaluate_columns(columns)

    def test_parse_peptidoform(self):
        test_cases = {
            "input": [
                ("_VGVGFGR_", 2),
                ("_MCK_", 2),
                ("_(ac)EEEIAALVIDNGSGMCK_", 2),
                ("_(gl)QYDADLEQILIQWITTQCRK_", 2),
                ("_LAM(ox)QEFMILPVGAANFR_", 2),
                ("_VGVN(de)GFGR_", 2),
                ("_(ac)EEEIAALVIDNGSGM(ox)CK_", 2),
                ("_(ac)SDKPDM(ox)AEIEK_", 2),
                ("_YYWGGHYSWDM(Ox)AK_", 2),
                ("_YYWGGHYSWDM(Oxidation (M))AK_", 2),
                ("_YYWGGHYM(ox)WDM(ox)AK_", 2),
                ("_(Acetyl (Protein N-term))ATGPM(ox)SFLK_", 2),
                ("_ACDE(Amidated (Peptide C-term))_", 2),
                ("_ACM(Ox)DE(Amidated (Peptide C-term))_", 2),
                # "_(Acetyl (Protein N-term))M(Ox)ACM(Ox)DEM(Ox)(Amidated (Peptide C-term))_",  # See levitsky/pyteomics/#77
            ],
            "expected_output": [
                "VGVGFGR/2",
                "MCK/2",
                "[ac]-EEEIAALVIDNGSGMCK/2",
                "[gl]-QYDADLEQILIQWITTQCRK/2",
                "LAM[ox]QEFMILPVGAANFR/2",
                "VGVN[de]GFGR/2",
                "[ac]-EEEIAALVIDNGSGM[ox]CK/2",
                "[ac]-SDKPDM[ox]AEIEK/2",
                "YYWGGHYSWDM[Ox]AK/2",
                "YYWGGHYSWDM[Oxidation (M)]AK/2",
                "YYWGGHYM[ox]WDM[ox]AK/2",
                "[Acetyl (Protein N-term)]-ATGPM[ox]SFLK/2",
                "ACDE-[Amidated (Peptide C-term)]/2",
                "ACM[Ox]DE-[Amidated (Peptide C-term)]/2",
                # "[Acetyl (Protein N-term)]-M[Ox]ACM[Ox]DEM[Ox]-[Amidated (Peptide C-term)]",
            ],
        }

        msms_reader = maxquant.MSMSReader("./tests/test_data/test_msms.txt")

        for test_in, expected_out in zip(test_cases["input"], test_cases["expected_output"]):
            output = msms_reader._parse_peptidoform(*test_in)
            assert output.proforma == expected_out


def test_issue_138_exact_row(caplog):
    from psm_utils.io import read_file

    with caplog.at_level("INFO", logger="psm_utils.io.maxquant"):
        [psm] = read_file("./tests/test_data/issue_138_msms.txt")
    assert psm.peptidoform.proforma == "ACVINGM[Oxidation]QLK/2"
    assert psm.peptidoform.theoretical_mass == pytest.approx(1091.5467367, abs=1e-4, rel=0)
    assert float(psm.metadata["Mass"]) - psm.peptidoform.theoretical_mass == pytest.approx(
        57.021464, abs=1e-4, rel=0
    )
    assert psm.spectrum_id == "18859"
    assert psm.precursor_mz == 575.29138
    assert psm.metadata["Sequence"] == "ACVINGMQLK"
    assert psm.metadata["Modifications"] == "Oxidation (M)"
    assert sum("modification mapping" in record.message for record in caplog.records) == 1


@pytest.fixture
def write_msms(tmp_path):
    import csv

    def write(*rows):
        path = tmp_path / "msms.txt"
        has_modifications = any(mods is not None for _, mods in rows)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(TEST_COL + (["Modifications"] if has_modifications else []))
            for scan, (sequence, mods) in enumerate(rows, start=1):
                row = ["run", scan, 2, sequence, "P1", 500, "", 20, 0.01, 100]
                writer.writerow(row + ([mods] if has_modifications else []))
        return path

    return write


@pytest.mark.parametrize(
    ("sequence", "modified", "expected", "delta"),
    [
        ("ACVINGMQLK", "_ACVINGM(Oxidation (M))QLK_", "ACVINGM[Oxidation]QLK", 15.994915),
        ("ACDEK", "_(Acetyl (Protein N-term))ACDEK_", "[Acetyl]-ACDEK", 42.010565),
        ("PEPTIDESK", "_PEPTIDES(Phospho (STY))K_", "PEPTIDES[Phospho]K", 79.966331),
        ("ACDK", "_AC(Carbamidomethyl (C))DK_", "AC[Carbamidomethyl]DK", 57.021464),
        ("NQDK", "_N(Deamidation (NQ))QDK_", "N[Deamidated]QDK", 0.984016),
        ("ACDEK", "_(Carbamyl (N-term))ACDEK_", "[Carbamyl]-ACDEK", 43.005814),
        ("PEPTIDEK", "_PEPTIDEK(GlyGly (K))_", "PEPTIDEK[GG]", 114.042927),
        ("ACDE", "_ACDE(Amidated (Peptide C-term))_", "ACDE-[Amidated]", -0.984016),
    ],
)
def test_descriptive_labels_resolve(write_msms, sequence, modified, expected, delta):
    from pyteomics import proforma

    from psm_utils import Peptidoform

    [psm] = maxquant.MSMSReader(write_msms((modified, None)))
    assert psm.peptidoform.proforma == f"{expected}/2"
    proforma.parse(psm.peptidoform.proforma)
    assert psm.peptidoform.sequence == sequence
    assert psm.peptidoform.theoretical_mass - Peptidoform(sequence).theoretical_mass == (
        pytest.approx(delta, abs=2e-6, rel=0)
    )


def test_short_labels_sites_repeats_counts_and_unmodified(write_msms):
    path = write_msms(
        (
            "_(ac)M(ox)AC(cm)N(de)M(ox)S(ph)_",
            "Acetyl (Protein N-term),2 Oxidation (M),Carbamidomethyl (C),Deamidation (NQ),Phospho (STY)",
        ),
        ("_(ca)AC(ca)DK_", "Carbamyl (N-term);Carbamidomethyl (C)"),
        ("_(gl)QPEPTIDEK(gl)_", "Gln->pyro-Glu,GlyGly (K)"),
        ("_(Gln->pyro-Glu)QPEPTIDE_", None),
        ("_(Glu->pyro-Glu)EPEPTIDE_", None),
        ("_ACDEK_", "Unmodified"),
    )
    assert [psm.peptidoform.proforma for psm in maxquant.MSMSReader(path)] == [
        "[Acetyl]-M[Oxidation]AC[Carbamidomethyl]N[Deamidated]M[Oxidation]S[Phospho]/2",
        "[Carbamyl]-AC[Carbamidomethyl]DK/2",
        "[Gln->pyro-Glu]-QPEPTIDEK[GG]/2",
        "[Gln->pyro-Glu]-QPEPTIDE/2",
        "[Glu->pyro-Glu]-EPEPTIDE/2",
        "ACDEK/2",
    ]


def test_same_token_at_two_sites_and_adjacent_terminal_tags(write_msms):
    path = write_msms(
        ("_(ac)ACDEK(ac)(Amidated (Peptide C-term))_", None),
    )
    [psm] = maxquant.MSMSReader(path)
    assert psm.peptidoform.proforma == "[Acetyl]-ACDEK[Acetyl]-[Amidated]/2"


def test_mapping_override_and_opt_out(write_msms, caplog):
    from psm_utils.io import read_file

    path = write_msms(("_(ac)AM(ox)C(in-house (C))DK_", None))
    overrides = {"ox": "UNIMOD:35", "in-house (C)": "+57.021464"}
    [raw] = read_file(path, modification_mapping={})
    assert raw.peptidoform.proforma == "[ac]-AM[ox]C[in-house (C)]DK/2"
    assert caplog.records == []
    [mapped] = read_file(path, modification_mapping=overrides)
    assert mapped.peptidoform.proforma == "[Acetyl]-AM[UNIMOD:35]C[+57.021464]DK/2"
    assert overrides == {"ox": "UNIMOD:35", "in-house (C)": "+57.021464"}


def test_unknown_labels_and_conflicting_short_code_warn_once(write_msms, caplog):
    path = write_msms(
        *[("_AM(ox)C(Novel.+ (tag (nested (C))))K_", "Oxidation variant (M)")] * 3,
        ("_AM(ox)CDEK_", "Oxidation (P)"),
        ("_PEPTIDEK(gl)_", "Glutaryl (K)"),
        ("_PEPTIDEK(gl)_", ""),
    )
    reader = maxquant.MSMSReader(path)
    for _ in range(2):
        assert [p.peptidoform.proforma for p in reader] == (
            ["AM[ox]C[Novel.+ (tag (nested (C)))]K/2"] * 3
            + ["AM[ox]CDEK/2", "PEPTIDEK[gl]/2", "PEPTIDEK[gl]/2"]
        )
    warnings = [
        r.message
        for r in caplog.records
        if r.name == "psm_utils.io.maxquant" and r.levelname == "WARNING"
    ]
    assert len(warnings) == 3
    assert any("'gl'" in msg for msg in warnings)
    assert any("'ox'" in msg and "modification_mapping" in msg for msg in warnings)
    assert any("'Novel.+ (tag (nested (C)))'" in msg for msg in warnings)


def test_already_valid_proforma_tags_are_preserved_without_warning(write_msms, caplog):
    path = write_msms(
        ("_(ac)AM(Unimod:35)C(+57.021464)K(U:1)_", None),
        ("_AK(Methyl)DEK_", None),
    )
    assert [p.peptidoform.proforma for p in maxquant.MSMSReader(path)] == [
        "[Acetyl]-AM[UNIMOD:35]C[+57.021464]K[UNIMOD:1]/2",
        "AK[Methyl]DEK/2",
    ]
    assert not [r for r in caplog.records if r.name == "psm_utils.io.maxquant"]
