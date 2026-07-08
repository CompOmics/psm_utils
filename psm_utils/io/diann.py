"""
Reader for PSM files from DIA-NN.

Reads the '.tsv' file as defined on the
`DIA-NN documentation page <https://github.com/vdemichev/DiaNN/tree/1.8.1?tab=readme-ov-file#main-output-reference>`_,
or the ``report.parquet`` main report written by DIA-NN 2.0 and later.

Notes
-----
- DIA-NN calculates q-values at both the run and library level. The run-level q-value is used as
  the PSM q-value.
- DIA-NN currently does not support C-terminal modifications in its searches.
- The classic ``.tsv`` main report does not contain precursor m/z values or decoy status; these
  are set to ``None``/``False``, respectively. The ``report.parquet`` main report does contain
  both, and they are used when available.
- The ``report.parquet`` main report does not contain a PSM-level score (no ``CScore`` column)
  or an MS2 scan number (no ``MS2.Scan`` column). The score is set to ``None`` and
  ``Precursor.Id`` is used as the spectrum identifier instead.

"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import]

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()

RESCORING_FEATURES: list[str] = [
    "RT",
    "Predicted.RT",
    "iRT",
    "Predicted.iRT",
    "Ms1.Profile.Corr",
    "Ms1.Area",
    "IM",
    "iIM",
    "Predicted.IM",
    "Predicted.iIM",
]


class DIANNTSVReader(ReaderBase):
    """Reader for DIA-NN TSV format."""

    def __init__(self, filename: str | Path, *args: Any, **kwargs: Any) -> None:
        """
        Reader for DIA-NN '.tsv' file.

        Parameters
        ----------
        filename : str or Path
            Path to PSM file.
        *args
            Additional positional arguments passed to the base class.
        **kwargs
            Additional keyword arguments passed to the base class.

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
        psm_dict: dict[str, str], filename: str | Path | None = None
    ) -> PSM:
        """Parse a single PSM from a DIA-NN PSM file."""
        rescoring_features: dict[str, Any] = {}
        for ft in RESCORING_FEATURES:
            try:
                rescoring_features[ft] = psm_dict[ft]
            except KeyError:
                continue

        return PSM(
            peptidoform=DIANNTSVReader._parse_peptidoform(
                psm_dict["Modified.Sequence"], psm_dict["Precursor.Charge"]
            ),
            spectrum_id=psm_dict["MS2.Scan"],
            run=psm_dict["Run"],
            is_decoy=False,
            qvalue=psm_dict["Q.Value"],
            pep=float(psm_dict["PEP"]),
            score=float(psm_dict["CScore"]),
            precursor_mz=None,  # Not returned by DIA-NN :(
            retention_time=float(psm_dict["RT"]),
            ion_mobility=float(psm_dict["IM"]),
            protein_list=psm_dict["Protein.Ids"].split(";"),
            source="diann",
            rank=None,
            provenance_data=({"diann_filename": str(filename)} if filename else {}),
            rescoring_features=rescoring_features,
            metadata={},
        )

    @staticmethod
    def _parse_peptidoform(peptide: str, charge: str | None) -> str:
        # Add charge
        if charge:
            peptide += f"/{int(float(charge))}"

        # Replace parentheses with square brackets and capitalize UniMod prefix
        pattern: str = r"\(UniMod:(\d+)\)"
        replacement: str = r"[UNIMOD:\1]"
        peptide = re.sub(pattern, replacement, peptide)

        # Add hyphen for N-terminal modifications
        # If [UNIMOD:n] occurs before the first amino acid, a hyphen is added before the first
        # amino acid
        if peptide[0] == "[":
            # Hyphen after the closing bracket
            peptide = peptide.replace("]", "]-", 1)

        # C-terminal modifications are currently not supported in DIA-NN

        return peptide

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame) -> PSMList:
        """Create a PSMList from a DIA-NN Pandas DataFrame."""
        records = cast(list[dict[str, str]], dataframe.to_dict(orient="records"))
        return PSMList(psm_list=[cls._get_peptide_spectrum_match(entry) for entry in records])


class DIANNParquetReader(ReaderBase):
    """Reader for the DIA-NN 2.0+ ``report.parquet`` main report."""

    def __init__(self, filename: str | Path, *args: Any, **kwargs: Any) -> None:
        """
        Reader for DIA-NN ``report.parquet`` file.

        Parameters
        ----------
        filename : str or Path
            Path to PSM file.
        *args
            Additional positional arguments passed to the base class.
        **kwargs
            Additional keyword arguments passed to the base class.

        """
        super().__init__(filename, *args, **kwargs)

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with pq.ParquetFile(self.filename) as pq_file:
            for batch in pq_file.iter_batches():
                for row in batch.to_pylist():
                    yield self._get_peptide_spectrum_match(row, self.filename)

    @staticmethod
    def _get_peptide_spectrum_match(
        psm_dict: dict[str, Any], filename: str | Path | None = None
    ) -> PSM:
        """Parse a single PSM from a DIA-NN ``report.parquet`` file."""
        rescoring_features: dict[str, Any] = {}
        for ft in RESCORING_FEATURES:
            try:
                rescoring_features[ft] = psm_dict[ft]
            except KeyError:
                continue

        return PSM(
            peptidoform=DIANNTSVReader._parse_peptidoform(
                psm_dict["Modified.Sequence"], str(psm_dict["Precursor.Charge"])
            ),
            spectrum_id=psm_dict["Precursor.Id"],
            run=psm_dict["Run"],
            is_decoy=bool(psm_dict["Decoy"]),
            qvalue=psm_dict["Q.Value"],
            pep=float(psm_dict["PEP"]),
            score=None,  # Not returned by DIA-NN in the `report.parquet` main report :(
            precursor_mz=float(psm_dict["Precursor.Mz"]),
            retention_time=float(psm_dict["RT"]),
            ion_mobility=float(psm_dict["IM"]),
            protein_list=psm_dict["Protein.Ids"].split(";"),
            source="diann",
            rank=None,
            provenance_data=({"diann_filename": str(filename)} if filename else {}),
            rescoring_features=rescoring_features,
            metadata={},
        )

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame) -> PSMList:
        """Create a PSMList from a DIA-NN ``report.parquet`` Pandas DataFrame."""
        records = cast(list[dict[str, Any]], dataframe.to_dict(orient="records"))
        return PSMList(psm_list=[cls._get_peptide_spectrum_match(entry) for entry in records])
