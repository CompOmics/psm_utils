"""Reader for PSM files from the AlphaDIA search engine."""

from __future__ import annotations

import csv
from abc import ABC
from collections.abc import Iterator
from pathlib import Path
from typing import Any, cast

import pandas as pd

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()

# TODO: check
RESCORING_FEATURES: list[str] = [
    "rt_observed",
    "mobility_observed",
    "mz_observed",
    "charge",
    "delta_rt",
]


class AlphaDIAReader(ReaderBase, ABC):
    """Reader for AlphaDIA TSV format."""

    def __init__(self, filename: str | Path, *args: Any, **kwargs: Any) -> None:
        """
        Reader for AlphaDIA ``precursor.tsv`` file.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments for parent class.
        **kwargs
            Additional keyword arguments for parent class.

        """
        super().__init__(filename, *args, **kwargs)

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")
            for row in reader:
                yield self._get_peptide_spectrum_match(row, self.filename)

    @staticmethod
    def _get_peptide_spectrum_match(
        psm_dict: dict[str, Any], filename: str | Path | None = None
    ) -> PSM:
        """Parse a single PSM from a AlphaDIA PSM file."""
        rescoring_features: dict[str, Any] = {}
        for ft in RESCORING_FEATURES:
            try:
                rescoring_features[ft] = psm_dict[ft]
            except KeyError:
                continue

        return PSM(
            peptidoform=AlphaDIAReader._parse_peptidoform(
                psm_dict["sequence"], psm_dict["mods"], psm_dict["mod_sites"], psm_dict["charge"]
            ),
            spectrum_id=psm_dict["frame_start"],  # TODO: needs to be checked
            run=psm_dict["run"],
            spectrum=psm_dict["frame_start"],  # TODO: needs to be checked
            is_decoy=bool(int(psm_dict["decoy"])),
            score=psm_dict["score"],
            qvalue=psm_dict["qval"],
            pep=psm_dict["proba"],
            precursor_mz=psm_dict["mz_observed"],
            retention_time=psm_dict["rt_observed"],
            ion_mobility=psm_dict["mobility_observed"],
            protein_list=psm_dict["proteins"].split(";"),
            rank=int(psm_dict["rank"]) + 1,  # AlphaDIA ranks are 0-based
            source="AlphaDIA",
            provenance_data=({"alphadia_filename": str(filename)} if filename else {}),
            metadata={},
            rescoring_features=rescoring_features,
        )

    @staticmethod
    def _parse_peptidoform(sequence: str, mods: str, mod_sites: str, charge: str | None) -> str:
        """Parse a peptidoform from a AlphaDIA PSM file."""
        # Parse modifications
        if mods:
            sequence_list: list[str] = [""] + list(sequence) + [""]  # N-term, sequence, C-term
            for mod, site_str in zip(mods.split(";"), mod_sites.split(";")):
                site: int = int(site_str)
                name: str = mod.split("@")[0]
                # N-terminal modification
                if site == 0:
                    sequence_list[0] = f"[{name}]-"
                # C-terminal modification
                elif site == -1:
                    sequence_list[-1] = f"-[{name}]"
                # Sequence modification
                else:
                    sequence_list[site] = f"{sequence_list[site]}[{name}]"
            sequence = "".join(sequence_list)

        # Add charge
        if charge:
            sequence += f"/{int(float(charge))}"

        return sequence

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame) -> PSMList:
        """Create a PSMList from a AlphaDIA Pandas DataFrame."""
        records = cast(list[dict[str, Any]], dataframe.to_dict(orient="records"))
        return PSMList(psm_list=[cls._get_peptide_spectrum_match(entry) for entry in records])
