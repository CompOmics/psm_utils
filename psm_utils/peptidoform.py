"""Peptidoform module for handling peptide sequences with modifications and charge states."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from typing import Literal, TypedDict, cast

import numpy as np
from pyteomics import mass, proforma  # type: ignore[import]
from pyteomics.proforma import ProForma, TagBase  # type: ignore[import]

from psm_utils.exceptions import PSMUtilsException
from psm_utils.utils import mass_to_mz


class Peptidoform:
    """Peptide sequence, modifications and charge state represented in ProForma notation."""

    def __init__(self, proforma_sequence: str | ProForma) -> None:
        """
        Peptide sequence, modifications and charge state represented in ProForma notation.

        Parameters
        ----------
        proforma_sequence
            Peptidoform sequence in ProForma v2 notation as :py:class:`str` or
            :py:class:`pyteomics.proforma.ProForma` object.

        Attributes
        ----------
        parsed_sequence : list
            List of tuples with residue and modifications for each location.
        properties : :py:class:`PeptidoformProperties`
            Dictionary with properties of the peptidoform, including N- and C-terminal
            modifications, unlocalized modifications, labile modifications, fixed
            modifications, and charge state.

        Examples
        --------
        >>> peptidoform = Peptidoform("ACDM[Oxidation]EK")
        >>> peptidoform.theoretical_mass
        711.2567622919099

        """
        self.parsed_sequence: list[tuple[str, list[TagBase] | None]]
        self.properties: PeptidoformProperties

        # Parse ProForma
        if isinstance(proforma_sequence, str):
            try:
                self.parsed_sequence, properties = proforma.parse(proforma_sequence)
                self.properties = cast(PeptidoformProperties, properties)
            except proforma.ProFormaError as e:
                raise PeptidoformException(
                    f"Could not parse ProForma sequence: {proforma_sequence}"
                ) from e
        elif isinstance(proforma_sequence, ProForma):
            self.parsed_sequence = proforma_sequence.sequence
            self.properties = proforma_sequence.properties
        else:
            raise TypeError(
                f"Expected ProForma sequence or ProForma object, got {type(proforma_sequence)}."
            )

        if self.properties["isotopes"]:
            raise NotImplementedError("Peptidoforms with isotopes are currently not supported.")

    def __repr__(self) -> str:
        """Return a string representation of the Peptidoform object."""
        return f"{self.__class__.__qualname__}('{self.proforma}')"

    def __str__(self) -> str:
        """Return the ProForma representation of the Peptidoform."""
        return self.proforma

    def __hash__(self) -> int:
        """Return a hash of the Peptidoform based on its ProForma representation."""
        return hash(self.proforma)

    def __eq__(self, __o: object) -> bool:
        """Check equality of Peptidoform with another object."""
        if isinstance(__o, str):
            return self.proforma == __o
        elif isinstance(__o, Peptidoform):  # type: ignore[return]
            return self.proforma == __o.proforma
        else:
            raise TypeError(f"Unsupported comparison type for Peptidoform: {type(__o)}")

    def __lt__(self, __o: object) -> bool:
        """Check if this Peptidoform is less than another object."""
        if isinstance(__o, str):
            return self.proforma < __o
        elif isinstance(__o, Peptidoform):
            return self.proforma < __o.proforma
        else:
            raise TypeError(f"Unsupported comparison type for Peptidoform: {type(__o)}")

    def __iter__(self) -> Iterable[tuple[str, None | list[TagBase]]]:
        """Return an iterator over the parsed sequence."""
        return self.parsed_sequence.__iter__()

    def __len__(self) -> int:
        """Return the length of the parsed sequence."""
        return self.parsed_sequence.__len__()

    def __getitem__(self, key: int) -> tuple[str, None | list[TagBase]]:
        """Get the item at the specified index from the parsed sequence."""
        return self.parsed_sequence.__getitem__(key)

    @property
    def proforma(self) -> str:
        """
        Peptidoform sequence in ProForma v2 notation.

        Examples
        --------
        >>> Peptidoform("AC[U:4]DEK/2").proforma
        'AC[UNIMOD:4]DEK/2'

        """
        return proforma.to_proforma(self.parsed_sequence, **self.properties)

    @property
    def sequence(self) -> str:
        """
        Stripped peptide sequence (modifications removed).

        Examples
        --------
        >>> Peptidoform("AC[U:4]DEK/2").sequence
        'ACDEK'

        """
        return "".join(pos[0] for pos in self.parsed_sequence)

    @property
    def modified_sequence(self) -> str:
        """
        Peptide sequence with modifications in ProForma format, but without charge state.

        Includes all modifications, including labile, unlocalized, and terminal modifications.

        Examples
        --------
        >>> Peptidoform("AC[U:4]DEK/2").modified_sequence
        'AC[U:4]DEK'

        """
        properties_without_charge = dict(self.properties).copy()
        properties_without_charge.pop("charge_state", None)
        return proforma.to_proforma(self.parsed_sequence, **properties_without_charge)

    @property
    def precursor_charge(self) -> int | None:
        """
        Returns the charge state as integer or :py:const:`None` if no charge assigned.

        Examples
        --------
        >>> Peptidoform("ACDEK/2").precursor_charge
        2

        >>> Peptidoform("ACDEK").precursor_charge
        None

        """
        try:
            return self.properties["charge_state"].charge
        except (AttributeError, KeyError):
            return None

    @property
    def is_modified(self) -> bool:
        """
        Whether or not the peptidoform carries any modification of any type.

        Includes N- and C-terminal, fixed, sequential, labile, and unlocalized
        modifications.

        """
        has_sequential = any(mods for _, mods in self.parsed_sequence)
        has_other = (
            bool(self.properties["n_term"])
            or bool(self.properties["c_term"])
            or bool(self.properties["unlocalized_modifications"])
            or bool(self.properties["labile_modifications"])
            or bool(self.properties["fixed_modifications"])
        )
        return has_sequential or has_other

    @property
    def sequential_composition(self) -> list[mass.Composition]:
        """
        Atomic compositions of both termini and each residue, including modifications.

        Includes N- and C-terminal, fixed, and sequential modifications. Does not
        include labile or unlocalized modifications.

        Examples
        --------
        >>> Peptidoform("ACDEK/2").sequential_composition
        [Composition({'H': 1}),
        Composition({'H': 5, 'C': 3, 'O': 1, 'N': 1}),
        Composition({'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1}),
        Composition({'H': 5, 'C': 4, 'O': 3, 'N': 1}),
        Composition({'H': 7, 'C': 5, 'O': 3, 'N': 1}),
        Composition({'H': 12, 'C': 6, 'N': 2, 'O': 1}),
        Composition({'H': 1, 'O': 1})]

        """
        # Get compositions for fixed modifications by amino acid
        fixed_rules = {}
        for rule in self.properties["fixed_modifications"]:
            if rule.targets is None:
                continue
            for aa in rule.targets:
                fixed_rules[aa] = rule.modification_tag.composition

        comp_list = []

        # N-terminus
        n_term = mass.Composition({"H": 1})
        if self.properties["n_term"]:
            for tag in self.properties["n_term"]:
                try:
                    n_term += tag.composition
                except (AttributeError, KeyError) as e:
                    raise ModificationException(
                        f"Cannot resolve composition for modification {tag.value}."
                    ) from e
        comp_list.append(n_term)

        # Sequence
        for aa, tags in self.parsed_sequence:
            # Amino acid
            try:
                position_comp = mass.std_aa_comp[aa].copy()
            except (AttributeError, KeyError) as e:
                # Allow X with modifications to specify gap with known composition
                if aa == "X":
                    if tags and all(hasattr(tag, "composition") for tag in tags):
                        position_comp = mass.Composition()
                    else:
                        raise AmbiguousResidueException(
                            "Cannot resolve composition for `X` without associated formula "
                            "modification."
                        ) from e
                else:
                    raise AmbiguousResidueException(
                        f"Cannot resolve composition for amino acid {aa}."
                    ) from e
            # Fixed modifications
            if aa in fixed_rules:
                position_comp += fixed_rules[aa]
            # Localized modifications
            if tags:
                for tag in tags:
                    tag = cast(proforma.ModificationBase, tag)
                    try:
                        position_comp += tag.composition
                    except (AttributeError, KeyError) as e:
                        raise ModificationException(
                            f"Cannot resolve composition for modification {tag.value}."
                        ) from e
            comp_list.append(position_comp)

        # C-terminus
        c_term = mass.Composition({"H": 1, "O": 1})
        if self.properties["c_term"]:
            for tag in self.properties["c_term"]:
                try:
                    c_term += tag.composition
                except (AttributeError, KeyError) as e:
                    raise ModificationException(
                        f"Cannot resolve composition for modification {tag.value}."
                    ) from e
        comp_list.append(c_term)

        return comp_list

    @property
    def composition(self) -> mass.Composition:
        """
        Atomic composition of the full peptidoform.

        Includes all modifications, also labile and unlocalized.

        Examples
        --------
        >>> Peptidoform("ACDEK/2").composition
        Composition({'H': 36, 'C': 21, 'O': 10, 'N': 6, 'S': 1})

        """
        comp = mass.Composition()
        for position_comp in self.sequential_composition:
            comp += position_comp
        for tag in self.properties["labile_modifications"]:
            try:
                comp += tag.composition
            except (AttributeError, KeyError) as e:
                raise ModificationException(
                    f"Cannot resolve composition for modification {tag.value}."
                ) from e
        for tag in self.properties["unlocalized_modifications"]:
            try:
                comp += tag.composition
            except (AttributeError, KeyError) as e:
                raise ModificationException(
                    f"Cannot resolve composition for modification {tag.value}."
                ) from e
        return comp

    @property
    def sequential_theoretical_mass(self) -> list[float]:
        """
        Monoisotopic mass of both termini and each (modified) residue.

        Includes N- and C-terminal, fixed, and sequential modifications. Does not
        include labile or unlocalized modifications.

        Examples
        --------
        >>> Peptidoform("ACDEK/2").sequential_theoretical_mass
        [1.00782503207,
        71.03711378471,
        103.00918478471,
        115.02694302383001,
        129.04259308796998,
        128.09496301399997,
        17.002739651629998]

        """
        fixed_rules = {}
        for rule in self.properties["fixed_modifications"]:
            if rule.targets is None:
                continue
            for aa in rule.targets:
                fixed_rules[aa] = rule.modification_tag.mass

        mass_list = []

        # N-terminus
        n_term = mass.Composition({"H": 1}).mass()
        if self.properties["n_term"]:
            for tag in self.properties["n_term"]:
                try:
                    n_term += tag.mass
                except (AttributeError, KeyError) as e:
                    raise ModificationException(
                        f"Cannot resolve mass for modification {tag.value}."
                    ) from e
        mass_list.append(n_term)

        # Sequence
        for aa, tags in self.parsed_sequence:
            # Amino acid
            try:
                position_mass = mass.std_aa_mass[aa]
            except (AttributeError, KeyError) as e:
                # Allow X with modifications to specify gap of unknown mass
                if aa == "X":
                    if tags and all(hasattr(tag, "mass") for tag in tags):
                        position_mass = 0.0
                    else:
                        raise AmbiguousResidueException(
                            "Cannot resolve mass for `X` without associated modification."
                        ) from e
                else:
                    raise AmbiguousResidueException(
                        f"Cannot resolve mass for amino acid {aa}."
                    ) from e
            # Fixed modifications
            if aa in fixed_rules:
                position_mass += fixed_rules[aa]
            # Localized modifications
            if tags:
                for tag in tags:
                    tag = cast(proforma.ModificationBase, tag)
                    try:
                        position_mass += tag.mass
                    except (AttributeError, KeyError) as e:
                        raise ModificationException(
                            f"Cannot resolve mass for modification {tag.value}."
                        ) from e
            mass_list.append(position_mass)

        # C-terminus
        c_term = mass.Composition({"H": 1, "O": 1}).mass()
        if self.properties["c_term"]:
            for tag in self.properties["c_term"]:
                try:
                    c_term += tag.mass
                except (AttributeError, KeyError) as e:
                    raise ModificationException(
                        f"Cannot resolve mass for modification {tag.value}."
                    ) from e
        mass_list.append(c_term)

        return mass_list

    @property
    def theoretical_mass(self) -> float:
        """
        Monoisotopic mass of the full uncharged peptidoform.

        Includes all modifications, also labile and unlocalized.

        Examples
        --------
        >>> Peptidoform("ACDEK/2").theoretical_mass
        564.22136237892

        """
        mass = sum(self.sequential_theoretical_mass)
        for tag in self.properties["labile_modifications"]:
            try:
                mass += tag.mass
            except (AttributeError, KeyError) as e:
                raise ModificationException(
                    f"Cannot resolve mass for modification {tag.value}."
                ) from e
        for tag in self.properties["unlocalized_modifications"]:
            try:
                mass += tag.mass
            except (AttributeError, KeyError) as e:
                raise ModificationException(
                    f"Cannot resolve mass for modification {tag.value}."
                ) from e
        return mass

    @property
    def theoretical_mz(self) -> float | None:
        """
        Monoisotopic mass-to-charge ratio of the full peptidoform.

        Includes all modifications, also labile and unlocalized.

        Examples
        --------
        >>> Peptidoform("ACDEK/2").theoretical_mz
        283.11850622153

        >>> Peptidoform("AC[+57.021464]DEK/2").theoretical_mz
        311.62923822153

        """
        if self.precursor_charge:
            return mass_to_mz(self.theoretical_mass, self.precursor_charge)
        else:
            return None

    def rename_modifications(self, mapping: dict[str, str]) -> None:
        """
        Apply mapping to rename modification tags.

        Parameters
        ----------
        mapping : dict[str, str]
            Mapping of ``old label`` → ``new label`` for each modification that
            requires renaming. Modification labels that are not in the mapping will not
            be renamed.

        See Also
        --------
        psm_utils.psm_list.PSMList.rename_modifications

        Examples
        --------
        >>> peptidoform = Peptidoform('[ac]-PEPTC[cmm]IDEK')
        >>> peptidoform.rename_modifications({
        ...     "ac": "Acetyl",
        ...     "cmm": "Carbamidomethyl"
        ... })
        >>> peptidoform.proforma
        '[Acetyl]-PEPTC[Carbamidomethyl]IDEK'

        """

        def _rename_modification_list(
            mods: list[TagBase] | None,
        ) -> list[TagBase] | None:
            if mods is None:
                return None

            new_mods = []
            for mod in mods:
                # Get value of the tag, formatted as string
                if isinstance(mod, proforma.MassModification):
                    mod_value = format_number_as_string(mod.value)
                else:
                    mod_value = mod.value

                # Rename modification if it is in the mapping
                if mod_value in mapping:
                    new_mods.append(proforma.process_tag_tokens(mapping[mod_value]))
                else:
                    new_mods.append(mod)

            return new_mods

        def _rename_modification_rule_list(
            rules: list[proforma.ModificationRule] | None,
        ) -> None:
            if rules is None:
                return None

            for rule in rules:
                # Rename modification tag if it is in the mapping
                if rule.modification_tag.value in mapping:
                    rule.modification_tag = proforma.process_tag_tokens(
                        mapping[rule.modification_tag.value]
                    )

        # Sequential modifications
        for i, (aa, mods) in enumerate(self.parsed_sequence):
            if mods:
                new_mods = _rename_modification_list(mods)
                self.parsed_sequence[i] = (aa, new_mods)

        # Non-sequence modifications
        for mod_type in (
            "n_term",
            "c_term",
            "unlocalized_modifications",
        ):
            self.properties[mod_type] = _rename_modification_list(self.properties[mod_type])  # type: ignore[assignment]

        # Modification rules
        _rename_modification_rule_list(self.properties["fixed_modifications"])

    def add_fixed_modifications(
        self, modification_rules: list[tuple[str, list[str]]] | dict[str, list[str]]
    ):
        """
        Add fixed modifications to peptidoform.

        Add modification rules for fixed modifications to peptidoform. These will be
        added in the "fixed modifications" notation, at the front of the ProForma
        sequence.

        See Also
        --------
        psm_utils.peptidoform.Peptidoform.apply_fixed_modifications

        Examples
        --------
        >>> peptidoform = Peptidoform("ATPEILTCNSIGCLK")
        >>> peptidoform.add_fixed_modifications([("Carbamidomethyl", ["C"])])
        >>> peptidoform.proforma
        '<[Carbamidomethyl]@C>ATPEILTCNSIGCLK'

        Notes
        -----
        While globally defined terminal modifications are not explicitly supported in ProForma v2,
        this function supports adding terminal modifications using the ``N-term`` and ``C-term``
        targets in place of an amino acid target. These global modifications are supported in the
        :py:meth:`psm_utils.peptidoform.Peptidoform.apply_fixed_modifications` method through a
        workaround. See https://github.com/HUPO-PSI/ProForma/issues/6 for discussions on the issue.

        """
        if isinstance(modification_rules, dict):
            modification_rules = list(modification_rules.items())

        parsed_modification_rules = [
            proforma.ModificationRule(proforma.process_tag_tokens(mod), targets)
            for mod, targets in modification_rules
        ]

        self.properties.setdefault("fixed_modifications", []).extend(parsed_modification_rules)

    def apply_fixed_modifications(self):
        """
        Apply ProForma fixed modifications as sequential modifications.

        Applies all modifications that are encoded as fixed in the ProForma notation
        (once at the beginning of the sequence) as modifications throughout the
        sequence at each affected amino acid residue.

        See Also
        --------
        psm_utils.peptidoform.Peptidoform.apply_fixed_modifications

        Examples
        --------
        >>> peptidoform = Peptidoform('<[Carbamidomethyl]@C>ATPEILTCNSIGCLK')
        >>> peptidoform.apply_fixed_modifications()
        >>> peptidoform.proforma
        'ATPEILTC[Carbamidomethyl]NSIGC[Carbamidomethyl]LK'

        """
        if self.properties["fixed_modifications"]:
            # Setup target_aa -> modification_list dictionary
            rule_dict = defaultdict(list)
            for rule in self.properties["fixed_modifications"]:
                if rule.targets is None:
                    continue
                for target_aa in rule.targets:
                    rule_dict[target_aa].append(rule.modification_tag)

            # Apply modifications to sequence
            for i, (aa, site_mods) in enumerate(self.parsed_sequence):
                if aa in rule_dict:
                    if site_mods:
                        self.parsed_sequence[i] = (aa, site_mods + rule_dict[aa])
                    else:
                        self.parsed_sequence[i] = (aa, rule_dict[aa])

            # Apply terminal modifications
            for term, term_name in [("n_term", "N-term"), ("c_term", "C-term")]:
                if term_name in rule_dict:
                    if self.properties[term]:
                        self.properties[term].extend(rule_dict[term_name])
                    else:
                        self.properties[term] = rule_dict[term_name]

            # Remove fixed modifications
            self.properties["fixed_modifications"] = []


class PeptidoformProperties(TypedDict):
    """Property items of a :py:class:`Peptidoform`."""

    n_term: list[proforma.ModificationBase] | None
    c_term: list[proforma.ModificationBase] | None
    unlocalized_modifications: list[proforma.ModificationBase]
    labile_modifications: list[proforma.ModificationBase]
    fixed_modifications: list[proforma.ModificationRule]
    charge_state: proforma.ChargeState
    isotopes: list[proforma.StableIsotope]


_ModificationsProperty = Literal[
    "n_term", "c_term", "unlocalized_modifications", "labile_modifications"
]

_ModificationRulesProperty = Literal["fixed_modifications"]


def format_number_as_string(num):
    """Format number as string for ProForma mass modifications."""
    # Using this method over `:+g` string formatting to avoid rounding and scientific notation
    num = float(num)
    plus = "+" if np.sign(num) == 1 else ""  # Add plus sign if positive
    num = str(num).rstrip("0").rstrip(".")  # Remove trailing zeros and decimal point
    return plus + num


class PeptidoformException(PSMUtilsException):
    """Error while handling :py:class:`Peptidoform`."""

    pass


class AmbiguousResidueException(PeptidoformException):
    """Error while handling ambiguous residue."""

    pass


class ModificationException(PeptidoformException):
    """Error while handling amino acid modification."""

    pass
