import pytest
from pyteomics import proforma

from psm_utils.peptidoform import Peptidoform, format_number_as_string


class TestPeptidoform:
    def test__len__(self):
        test_cases = [
            ("ACDEFGHIK", 9),
            ("[ac]-AC[cm]DEFGHIK", 9),
            ("[ac]-AC[Carbamidomethyl]DEFGHIK", 9),
            ("[Acetyl]-AC[cm]DEFGK", 7),
            ("<[cm]@C>[Acetyl]-ACDK", 4),
            ("<[Carbamidomethyl]@C>[ac]-ACDEFGHIK", 9),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert len(peptidoform) == expected_out

    def test__eq__(self):
        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK", True),
            ("ACDEFGHIK", "ACDEFGHI", False),
            ("ACDEFGHIK/2", "ACDEFGHIK/2", True),
            ("ACDEFGHIK/2", "ACDEFGHIK/3", False),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGHIK", True),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGH", False),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGH", False),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGH", False),
        ]

        for test_case_in_1, test_case_in_2, expected_out in test_cases:
            assert (Peptidoform(test_case_in_1) == test_case_in_2) == expected_out
            assert (Peptidoform(test_case_in_1) == Peptidoform(test_case_in_2)) == expected_out

    with pytest.raises(TypeError):
        Peptidoform("ACDEFGHIK") == 1

    def test__getitem__(self):
        test_cases = [
            ("ACDEFGHIK", 0, ("A", None)),
            ("ACDEFGHIK", 8, ("K", None)),
            ("[ac]-AC[cm]DEFGHIK", 0, ("A", None)),
            ("[ac]-AC[cm]DEFGHIK", 1, ("C", [proforma.GenericModification("cm")])),
            ("[ac]-AC[cm]DEFGHIK", 8, ("K", None)),
        ]

        for test_case_in, index, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert peptidoform[index] == expected_out

    def test__iter__(self):
        for aa, mods in Peptidoform("ACDEM[U:35]K"):
            assert isinstance(aa, str)
            if mods is not None:
                assert isinstance(mods, list)
                for mod in mods:
                    assert isinstance(mod, proforma.TagBase)

    def test_sequence(self):
        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "ACDEFGHIK"),
            ("[ac]-AC[Carbamidomethyl]DEFGHIK", "ACDEFGHIK"),
            ("[Acetyl]-AC[cm]DEFGK", "ACDEFGK"),
            ("<[cm]@C>[Acetyl]-ACDK", "ACDK"),
            ("<[Carbamidomethyl]@C>[ac]-ACDEFGHIK", "ACDEFGHIK"),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert peptidoform.sequence == expected_out

    def test_modified_sequence(self):
        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK"),
            ("ACDEFGHIK/3", "ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK/3", "[ac]-AC[cm]DEFGHIK"),
            ("<[cm]@C>[Acetyl]-ACDK/3", "<[cm]@C>[Acetyl]-ACDK"),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert peptidoform.modified_sequence == expected_out

    def test_precursor_charge(self):
        test_cases = [
            ("ACDEFGHIK", None),
            ("ACDEFGHIK/2", 2),
            ("ACDEFGHIK/3", 3),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert peptidoform.precursor_charge == expected_out

    def test_rename_modifications(self):
        label_mapping = {
            "ac": "Acetyl",
            "cm": "Carbamidomethyl",
            "+57.021": "Carbamidomethyl",
            "-18.010565": "Glu->pyro-Glu",
        }

        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("[ac]-AC[Carbamidomethyl]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("[Acetyl]-AC[cm]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("<[cm]@C>[Acetyl]-ACDEFGHIK", "<[Carbamidomethyl]@C>[Acetyl]-ACDEFGHIK"),
            ("<[Carbamidomethyl]@C>[ac]-ACDEFGHIK", "<[Carbamidomethyl]@C>[Acetyl]-ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("AC[+57.021]DEFGHIK", "AC[Carbamidomethyl]DEFGHIK"),
            ("E[-18.010565]DEK", "E[Glu->pyro-Glu]DEK"),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            peptidoform.rename_modifications(label_mapping)
            assert peptidoform.proforma == expected_out

    def test_add_apply_fixed_modifications(self):
        test_cases = [
            ("ACDEK", [("Cmm", ["C"])], "AC[Cmm]DEK"),
            ("AC[Cmm]DEK", [("SecondMod", ["C"])], "AC[Cmm][SecondMod]DEK"),
            ("ACDEK", [("TMT6plex", ["K", "N-term"])], "[TMT6plex]-ACDEK[TMT6plex]"),
            ("ACDEK-[CT]", [("TMT6plex", ["K", "N-term"])], "[TMT6plex]-ACDEK[TMT6plex]-[CT]"),
        ]

        for test_case_in, fixed_modifications, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            peptidoform.add_fixed_modifications(fixed_modifications)
            peptidoform.apply_fixed_modifications()
            assert peptidoform.proforma == expected_out

    def test_sequential_theoretical_mass(self):
        """Test sequential theoretical mass calculation."""
        test_cases = [
            # Simple peptide: (proforma_str, number_of_residues)
            ("ACDEK", 5),  # N-term, A, C, D, E, K, C-term = 7 total
            # Peptide with modifications
            ("[Acetyl]-ACDEK", 5),
            ("AC[Carbamidomethyl]DEK", 5),
            # Peptide with X and mass modification (gap of known mass)
            ("ACX[+100.5]DEK", 6),  # A, C, X, D, E, K
            ("X[+50.0]ACDE", 5),  # X, A, C, D, E
            # Multiple X residues with mass modifications
            ("X[+100.0]ACX[+200.0]DE", 6),  # X, A, C, X, D, E
        ]

        for proforma_str, num_residues in test_cases:
            peptidoform = Peptidoform(proforma_str)
            seq_mass = peptidoform.sequential_theoretical_mass

            # Check that we get the right number of elements (N-term + residues + C-term)
            expected_length = num_residues + 2  # +2 for N-term and C-term
            assert len(seq_mass) == expected_length, (
                f"Failed for {proforma_str}: expected {expected_length}, got {len(seq_mass)}"
            )

            # Check that all values are floats
            assert all(isinstance(m, float) for m in seq_mass), f"Failed for {proforma_str}"

            # Check that sum matches theoretical mass (excluding charge)
            total_mass = sum(seq_mass)
            expected_total = peptidoform.theoretical_mass
            assert abs(total_mass - expected_total) < 1e-6, (
                f"Failed for {proforma_str}: {total_mass} != {expected_total}"
            )

    def test_sequential_theoretical_mass_with_x_gap(self):
        """Test sequential theoretical mass with X representing a gap of known mass."""
        # X[+100.5] should contribute 100.5 to the mass
        peptidoform = Peptidoform("ACX[+100.5]DE")
        seq_mass = peptidoform.sequential_theoretical_mass

        # seq_mass should be: [N-term, A, C, X+100.5, D, E, C-term]
        assert len(seq_mass) == 7

        # The X residue (index 3) should have mass 0.0 + 100.5 = 100.5
        x_mass = seq_mass[3]
        assert abs(x_mass - 100.5) < 1e-6, f"Expected 100.5, got {x_mass}"

    def test_sequential_theoretical_mass_with_x_no_modification_fails(self):
        """Test that X without modification fails for sequential_theoretical_mass."""
        from psm_utils.peptidoform import AmbiguousResidueException

        # X without any modification should fail for mass calculation
        peptidoform = Peptidoform("ACXDE")

        with pytest.raises(
            AmbiguousResidueException,
            match="Cannot resolve mass for `X` without associated modification",
        ):
            _ = peptidoform.sequential_theoretical_mass

    def test_sequential_composition(self):
        """Test sequential composition calculation."""
        from pyteomics import mass

        test_cases = [
            # Simple peptide: (proforma_str, number_of_residues)
            ("ACDEK", 5),  # N-term, A, C, D, E, K, C-term = 7 total
            # Peptide with modifications
            ("[Acetyl]-ACDEK", 5),
            ("AC[Carbamidomethyl]DEK", 5),
            # Peptide with terminal modifications
            ("[Acetyl]-ACDEK-[Amidated]", 5),
        ]

        for proforma_str, num_residues in test_cases:
            peptidoform = Peptidoform(proforma_str)
            seq_comp = peptidoform.sequential_composition

            # Check that we get the right number of elements (N-term + residues + C-term)
            expected_length = num_residues + 2  # +2 for N-term and C-term
            assert len(seq_comp) == expected_length, (
                f"Failed for {proforma_str}: expected {expected_length}, got {len(seq_comp)}"
            )

            # Check that all values are Composition objects
            assert all(isinstance(c, mass.Composition) for c in seq_comp), (
                f"Failed for {proforma_str}"
            )

            # Check that sum matches full composition
            total_comp = mass.Composition()
            for comp in seq_comp:
                total_comp += comp
            assert total_comp == peptidoform.composition, f"Failed for {proforma_str}"

    def test_sequential_composition_with_x_gap(self):
        """Test sequential composition with X representing a gap of unknown composition."""
        from pyteomics import mass

        # X with formula modification should allow empty base composition
        peptidoform = Peptidoform("ACX[Formula:C6H12O6]DE")
        seq_comp = peptidoform.sequential_composition

        # seq_comp should be: [N-term, A, C, X+composition, D, E, C-term]
        assert len(seq_comp) == 7

        # The X residue (index 3) should have composition C6H12O6
        x_comp = seq_comp[3]
        expected_comp = mass.Composition({"C": 6, "H": 12, "O": 6})
        assert x_comp == expected_comp, f"Expected {expected_comp}, got {x_comp}"

    def test_sequential_composition_with_x_mass_only_fails(self):
        """Test that X with only mass modification fails for sequential_composition."""
        from psm_utils.peptidoform import AmbiguousResidueException

        # X with only mass modification should fail for composition calculation
        peptidoform = Peptidoform("ACX[+100.5]DE")

        with pytest.raises(
            AmbiguousResidueException,
            match="Cannot resolve composition for `X` without associated formula modification",
        ):
            _ = peptidoform.sequential_composition


def test_format_number_as_string():
    """Test format_number_as_string function."""
    test_cases = [
        (1212.12, "+1212.12"),
        (-1212.12, "-1212.12"),
        (0.1, "+0.1"),
        (-0.1, "-0.1"),
        (1212.000, "+1212"),
        (1212.1200, "+1212.12"),
    ]

    for test_case_in, expected_out in test_cases:
        assert format_number_as_string(test_case_in) == expected_out
