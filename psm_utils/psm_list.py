from __future__ import annotations

import re
from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd
from pydantic import BaseModel
from pyteomics import auxiliary, proforma
from rich.pretty import pretty_repr

from psm_utils.psm import NUMPY_DTYPES, PSM


class PSMList(BaseModel):
    """Data class representing a list of PSMs."""

    psm_list: List[PSM]

    def __init__(__pydantic_self__, **data) -> None:
        """
        Data class representing a list of PSMs, with some useful functionality.

        Parameters
        ----------
        psm_list : list[PSM]
            List of PSM instances.

        Examples
        --------
        Initiate a :py:class:`PSMList` from a list of PSM objects:

        >>> psm_list = PSMList(psm_list=[
        ...     PSM(peptidoform="ACDK", spectrum_id=1, score=140.2, retention_time=600.2),
        ...     PSM(peptidoform="CDEFR", spectrum_id=2, score=132.9, retention_time=1225.4),
        ...     PSM(peptidoform="DEM[Oxidation]K", spectrum_id=3, score=55.7, retention_time=3389.1),
        ... ])

        :py:class:`PSMList` directly supports iteration:

        >>> for psm in psm_list:
        ...     print(psm.peptidoform.score)
        140.2
        132.9
        55.7

        :py:class:`PSM` properties can be accessed as a single Numpy array:

        >>> psm_list["score"]
        array([140.2, 132.9, 55.7], dtype=object)

        :py:class:`PSMList` supports indexing and slicing:

        >>> psm_list_subset = psm_list[0:2]
        >>> psm_list_subset["score"]
        array([140.2, 132.9], dtype=object)

        >>> psm_list_subset = psm_list[0, 2]
        >>> psm_list_subset["score"]
        array([140.2, 55.7], dtype=object)

        For more advanced and efficient vectorized access, converting the
        :py:class:`PSMList` to a Pandas DataFrame is highly recommended:

        >>> psm_df = psm_list.to_dataframe()
        >>> psm_df[(psm_df["retention_time"] < 2000) & (psm_df["score"] > 10)]
          peptidoform  spectrum_id   run collection spectrum is_decoy  score qvalue   pep precursor_mz  retention_time protein_list  rank source provenance_data metadata rescoring_features
        0        ACDK            1  None       None     None     None  140.2   None  None         None           600.0         None  None   None            None     None               None
        1       CDEFR            2  None       None     None     None  132.9   None  None         None          1225.0         None  None   None            None     None               None

        """
        super().__init__(**data)

    def __rich_repr__(self):
        yield "psm_list", self.psm_list

    def __repr__(self):
        return pretty_repr(self, max_length=5)

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        return PSMList(psm_list=self.psm_list + other.psm_list)

    def __iter__(self) -> Iterable[PSM]:
        return self.psm_list.__iter__()

    def __len__(self) -> int:
        return self.psm_list.__len__()

    def __getitem__(self, item) -> PSM | list[PSM]:
        if isinstance(item, (int, np.integer)):
            # Return single PSM by index
            return self.psm_list[item]
        elif isinstance(item, slice):
            # Return new PSMList from slice
            return PSMList(psm_list=self.psm_list[item])
        elif isinstance(item, str):
            # Return PSM property as array across full PSMList
            try:
                return np.fromiter(
                    (psm[item] for psm in self.psm_list), dtype=NUMPY_DTYPES[item], count=len(self)
                )
            except TypeError:
                return np.fromiter(
                    (psm[item] for psm in self.psm_list), dtype=object, count=len(self)
                )
        elif _is_iterable_of_bools(item):
            # Return new PSMList with items that were True
            return PSMList(psm_list=[self.psm_list[i] for i in np.flatnonzero(item)])
        elif _is_iterable_of_ints(item):
            # Return new PSMList with selection of PSMs by list indices
            return PSMList(psm_list=[self.psm_list[i] for i in item])
        elif _is_iterable_of_strings(item):
            # Return 2D numpy array of multiple PSM properties along full PSMList
            return np.array([self[inner_item] for inner_item in item]).transpose()
        else:
            raise TypeError(f"Unsupported indexing type: {type(item)}")

    def __setitem__(self, item, values: Sequence) -> None:
        if not len(values) == len(self):
            raise ValueError(f"Expected value with same length as PSMList: {len(self)}")
        for value, psm in zip(values, self):
            psm[item] = value

    @property
    def collections(self) -> list:
        """List of collections in :py:class:`PSMList`."""
        if (self["collection"] != None).any():  # noqa: E711
            return list(np.unique(self["collection"]))
        else:
            return [None]

    @property
    def runs(self) -> list:
        """List of runs in :py:class:`PSMList`."""
        if (self["run"] != None).any():  # noqa: E711
            return list(np.unique(self["run"]))
        else:
            return [None]

    def append(self, psm: PSM) -> None:
        """Append PSM to :py:class:`PSMList`."""
        self.psm_list.append(psm)

    def extend(self, psm_list: PSMList) -> None:
        """Extend :py:class:`PSMList` with another :py:class:`PSMList`."""
        self.psm_list.extend(psm_list)

    def get_psm_dict(self):
        """Get nested dictionary of PSMs by collection, run, and spectrum_id."""
        psm_dict = {}
        for psm in self:
            try:
                psm_dict[psm.collection][psm.run][psm.spectrum_id].append(psm)
            except KeyError:
                try:
                    psm_dict[psm.collection][psm.run][psm.spectrum_id] = [psm]
                except KeyError:
                    try:
                        psm_dict[psm.collection][psm.run] = {psm.spectrum_id: [psm]}
                    except KeyError:
                        psm_dict[psm.collection] = {psm.run: {psm.spectrum_id: [psm]}}
        return psm_dict

    def set_ranks(self, lower_score_better: bool = False):
        """Set identification ranks for all PSMs in :py:class:`PSMList`."""
        columns = ["collection", "run", "spectrum_id", "score"]
        self["rank"] = (
            pd.DataFrame(self[columns], columns=columns)
            .sort_values("score", ascending=lower_score_better)
            .fillna(0)  # groupby does not play well with None values
            .groupby(["collection", "run", "spectrum_id"])
            .cumcount()
            .sort_index()
            + 1  # 1-based counting
        )

    def get_rank1_psms(self, *args, **kwargs) -> PSMList:
        """
        Return new :py:class:`PSMList` with only first-rank PSMs.

        First runs :py:meth:`~set_ranks` with ``*args`` and ``**kwargs`` if if any PSM
        has no rank yet.
        """
        if None in self["rank"]:
            self.set_ranks(*args, **kwargs)
        return self[self["rank"] == 1]

    def find_decoys(self, decoy_pattern: str) -> None:
        """
        Use regular expression pattern to find decoy PSMs by protein name(s).

        This method allows a regular expression pattern to be applied on
        :py:obj:`~psm_utils.psm.PSM`
        :py:attr:`~psm_utils.psm.PSM.protein_list` items to set the
        :py:attr:`~psm_utils.psm.PSM.is_decoy` attribute.
        Decoy protein entries are commonly marked with a prefix or suffix, e.g.
        ``DECOY_``, or ``_REVERSED``. If ``decoy_pattern`` matches to a substring of all
        entries in :py:attr:`protein_list`, the PSM is interpreted as a decoy. Existing
        :py:attr:`is_decoy` entries are overwritten.

        Parameters
        ----------
        decoy_pattern: str
            Regular expression pattern to match decoy protein entries.

        Examples
        --------
        >>> psm_list.find_decoys(r"^DECOY_")

        """
        decoy_pattern = re.compile(decoy_pattern)
        for psm in self:
            psm.is_decoy = all([decoy_pattern.search(p) is not None for p in psm.protein_list])

    def calculate_qvalues(self, reverse: bool = True, **kwargs) -> None:
        """
        Calculate q-values using the target-decoy approach.

        Q-values are calculated for all PSMs from the target and decoy scores. This
        requires that all PSMs have a :py:attr:`score` and a target/decoy state
        (:py:attr:`is_decoy`) assigned. Any existing q-values will be overwritten.

        Parameters
        ----------
        reverse: boolean, optional
            If True (default), a higher score value indicates a better PSM.
        **kwargs : dict, optional
            Additional arguments to be passed to
            `pyteomics.auxiliary.target_decoy.qvalues <https://pyteomics.readthedocs.io/en/latest/api/auxiliary.html#pyteomics.auxiliary.target_decoy.qvalues>`_.

        """
        for key in ["score", "is_decoy"]:
            if (self[key] == None).any():  # noqa: E711 (self[key] is a Numpy array)
                raise ValueError(
                    f"Cannot calculate q-values if not all PSMs have `{key}` assigned."
                )

        qvalues = auxiliary.qvalues(
            self,
            key="score",
            is_decoy="is_decoy",
            formula=1,
            correction=1,
            remove_decoy=False,
            reverse=reverse,
            full_output=True,
            **kwargs,
        )
        for score, is_decoy, qvalue, psm in qvalues:
            psm.qvalue = qvalue

    def rename_modifications(self, mapping: dict[str, str]) -> None:
        """
        Apply mapping to rename modification tags for all PSMs.

        Applies :py:meth:`psm_utils.peptidoform.Peptidoform.rename_modifications` on
        all PSM peptidoforms in the :py:class:`PSMList`.

        Parameters
        ----------
        mapping : dict[str, str]
            Mapping of ``old label`` → ``new label`` for each modification that
            requires renaming. Modification labels that are not in the mapping will not
            be renamed.

        See also
        --------
        psm_utils.peptidoform.Peptidoform.rename_modifications

        """
        for psm in self.psm_list:
            psm.peptidoform.rename_modifications(mapping)

    def add_fixed_modifications(
        self, modification_rules: list[tuple[str, list[str]]] | dict[str, list[str]]
    ):
        """
        Add fixed modifications to all PSM peptidoforms in :py:class:`PSMList`.

        Add modification rules for fixed modifications to peptidoform. These will be
        added in the "fixed modifications" notation, at the front of the ProForma
        sequence.

        See also
        --------
        psm_utils.peptidoform.Peptidoform.add_fixed_modifications

        Examples
        --------
        >>> psm_list.add_fixed_modifications([("Carbamidomethyl", ["C"])])

        >>> psm_list.add_fixed_modifications({"Carbamidomethyl": ["C"]})

        """
        if isinstance(modification_rules, dict):
            modification_rules = modification_rules.items()
        modification_rules = [
            proforma.ModificationRule(proforma.process_tag_tokens(mod), targets)
            for mod, targets in modification_rules
        ]
        for psm in self.psm_list:
            if psm.peptidoform.properties["fixed_modifications"]:
                psm.peptidoform.properties["fixed_modifications"].extend(modification_rules)
            else:
                psm.peptidoform.properties["fixed_modifications"] = modification_rules

    def apply_fixed_modifications(self):
        """
        Apply ProForma fixed modifications as sequential modifications.

        Applies :py:meth:`psm_utils.peptidoform.Peptidoform.apply_fixed_modifications`
        on all PSM peptidoforms in the :py:class:`PSMList`.

        See also
        --------
        psm_utils.peptidoform.Peptidoform.apply_fixed_modifications

        Examples
        --------
        >>> psm_list.apply_fixed_modifications()

        """
        for psm in self.psm_list:
            psm.peptidoform.apply_fixed_modifications()

    def to_dataframe(self) -> pd.DataFrame:
        """Convert :py:class:`PSMList` to :py:class:`pandas.DataFrame`."""
        return pd.DataFrame.from_records([psm.__dict__ for psm in self])


def _is_iterable_of_bools(obj):
    try:
        if all(isinstance(x, (bool, np.bool_)) for x in obj):
            return True
        else:
            return False
    except (TypeError, ValueError):
        return False


def _is_iterable_of_ints(obj):
    try:
        if not all(isinstance(x, (int, np.integer)) for x in obj):
            return False
        else:
            return True
    except (TypeError, ValueError):
        return False


def _is_iterable_of_strings(obj):
    try:
        if not all(isinstance(x, str) for x in obj):
            return False
        else:
            return True
    except (TypeError, ValueError):
        return False
