"""Reader for ProteoScape Parquet files."""

import logging
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import]

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import format_number_as_string
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)

DECOY_PATTERN = re.compile(r"^Reverse_")


class ProteoScapeReader(ReaderBase):
    """Reader for ProteoScape Parquet files."""

    def __init__(
        self,
        filename: str | Path,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """
        Reader for ProteoScape Parquet files.

        Parameters
        ----------
        filename
            Path to ProteoScape Parquet file.
        *args
            Additional positional arguments passed to the base class.
        **kwargs
            Additional keyword arguments passed to the base class.

        """
        super().__init__(filename, *args, **kwargs)

    def __len__(self) -> int:
        """Return number of PSMs in file."""
        return pq.read_metadata(self.filename).num_rows

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with pq.ParquetFile(self.filename) as reader:
            for batch in reader.iter_batches():
                for row in batch.to_pylist():
                    try:
                        yield _parse_entry(row)
                    except Exception as e:
                        raise PSMUtilsIOException(f"Error while parsing row {row}:\n{e}") from e

    @staticmethod
    def from_dataframe(dataframe: pd.DataFrame) -> PSMList:
        """Create a PSMList from a ProteoScape Pandas DataFrame."""
        return PSMList(
            psm_list=[_parse_entry(entry) for entry in dataframe.to_dict(orient="records")]  # type: ignore[arg-type]
        )


def _parse_peptidoform(
    stripped_peptide: str,
    ptms: np.ndarray[Any, Any],
    ptm_locations: np.ndarray[Any, Any],
    precursor_charge: int,
) -> str:
    """Parse peptide sequence and modifications to ProForma."""
    peptidoform = list(stripped_peptide)
    n_term = ""
    c_term = ""
    for ptm, ptm_location in zip(ptms, ptm_locations, strict=True):
        ptm_str = format_number_as_string(ptm)
        if ptm_location == -1:
            n_term = f"[{ptm_str}]-"
        elif ptm_location == len(peptidoform):
            c_term = f"-[{ptm_str}]"
        else:
            peptidoform[ptm_location] = f"{peptidoform[ptm_location]}[{ptm_str}]"
    return f"{n_term}{''.join(peptidoform)}{c_term}/{precursor_charge}"


def _parse_entry(entry: dict[str, Any]) -> PSM:
    """Parse a single entry from ProteoScape Parquet file to PSM object."""
    return PSM(
        peptidoform=_parse_peptidoform(
            entry["stripped_peptide"],
            entry["ptms"],
            entry["ptm_locations"],
            entry["precursor_charge"],
        ),
        spectrum_id=entry["ms2_id"],
        run=entry.get("run"),
        is_decoy=all(DECOY_PATTERN.match(p) for p in entry["locus_name"]),
        score=entry["x_corr_score"],
        precursor_mz=entry["precursor_mz"],
        retention_time=entry["rt"],
        ion_mobility=entry["ook0"],
        protein_list=list(entry["locus_name"]),
        rank=entry["rank"],
        source="ProteoScape",
        provenance_data={
            "candidate_id": str(entry["candidate_id"]),
            "ms2_id": str(entry["ms2_id"]),
            "parent_id": str(entry["parent_id"]),
        },
        metadata={
            "leading_aa": str(entry["leading_aa"]),
            "trailing_aa": str(entry["trailing_aa"]),
            "corrected_ook0": str(entry["corrected_ook0"]),
        },
        rescoring_features={
            "tims_score": float(entry["tims_score"]),
            "x_corr_score": float(entry["x_corr_score"]),
            "delta_cn_score": float(entry["delta_cn_score"]),
            "ppm_error": float(entry["ppm_error"]),
            "number_matched_ions": float(entry["number_matched_ions"]),
            "number_expected_ions": float(entry["number_expected_ions"]),
            "ion_proportion": float(entry["ion_proportion"]),
            "spectrum_total_ion_intensity": float(entry["spectrum_total_ion_intensity"]),
        },
    )
