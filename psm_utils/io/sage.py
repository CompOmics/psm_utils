"""
Reader for PSM files from the Sage search engine.

Reads the ``results.sage.tsv`` file as defined on the
`Sage documentation page <https://github.com/lazear/sage/blob/v0.14.7/DOCS.md#interpreting-sage-output>`_.

"""

from __future__ import annotations

import csv
from abc import ABC, abstractmethod
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import]
from pyteomics import mass  # type: ignore[import]

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()


class _SageReaderBase(ReaderBase, ABC):
    def __init__(
        self,
        filename: str | Path,
        *args: Any,
        score_column: str = "sage_discriminant_score",
        **kwargs: Any,
    ) -> None:
        """
        Reader for Sage results file.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments passed to parent class.
        score_column
            Name of the column that holds the primary PSM score. Default is
            ``sage_discriminant_score``, ``hyperscore`` could also be used.
        **kwargs
            Additional keyword arguments passed to parent class.

        """
        super().__init__(filename, *args, **kwargs)
        self.score_column = score_column

    @abstractmethod
    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        raise NotImplementedError("Use `SageTSVReader` or `SageParquetReader` instead.")

    def _get_peptide_spectrum_match(self, psm_dict: dict[str, Any]) -> PSM:
        """Parse a single PSM from a sage PSM file."""
        rescoring_features: dict[str, Any] = {}
        for ft in RESCORING_FEATURES:
            try:
                rescoring_features[ft] = psm_dict[ft]
            except KeyError:
                continue

        ion_mobility_features = self._extract_ion_mobility_features(psm_dict)
        rescoring_features.update(ion_mobility_features)

        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["peptide"],
                psm_dict["charge"],
            ),
            spectrum_id=psm_dict["scannr"],
            run=Path(psm_dict["filename"]).stem,
            is_decoy=psm_dict["is_decoy"],
            qvalue=psm_dict["spectrum_q"],
            score=float(psm_dict[self.score_column]),
            precursor_mz=self._parse_precursor_mz(psm_dict["expmass"], psm_dict["charge"]),
            retention_time=float(psm_dict["rt"]),
            ion_mobility=rescoring_features.get("ion_mobility"),
            protein_list=psm_dict["proteins"].split(";"),
            source="sage",
            rank=int(float(psm_dict["rank"])),
            provenance_data={"sage_filename": self.filename.as_posix()},
            rescoring_features=rescoring_features,
            metadata={},
        )

    @staticmethod
    def _parse_peptidoform(peptide: str, charge: str | None) -> str:
        """Parse peptide sequence and charge to peptidoform string."""
        if charge:
            peptide += f"/{int(float(charge))}"
        return peptide

    @staticmethod
    def _parse_precursor_mz(expmass: str, charge: str | None) -> float | None:
        """Parse experimental mass and charge to precursor m/z."""
        if charge:
            charge_val = float(charge)
            expmass_val = float(expmass)
            return (expmass_val + (mass.nist_mass["H"][1][0] * charge_val)) / charge_val
        return None

    @staticmethod
    def _extract_ion_mobility_features(psm_dict: dict[str, Any]) -> dict[str, float]:
        """Extract ion mobility features from the PSM dictionary if present and non-zero."""
        try:
            ion_mob = float(psm_dict["ion_mobility"])
            if ion_mob:
                return {
                    "ion_mobility": ion_mob,
                    "predicted_mobility": float(psm_dict["predicted_mobility"]),
                    "delta_mobility": float(psm_dict["delta_mobility"]),
                }
        except (KeyError, ValueError):
            pass
        return {}

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame) -> PSMList:
        """Create a PSMList from a Sage Pandas DataFrame."""
        return PSMList(
            psm_list=[
                cls._get_peptide_spectrum_match(cls(""), entry)  # type: ignore[arg-type]
                for entry in dataframe.to_dict(orient="records")
            ]
        )


class SageTSVReader(_SageReaderBase):
    """Reader for Sage TSV results files."""

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as open_file:
            reader = csv.DictReader(open_file, delimiter="\t")
            for row in reader:
                row["is_decoy"] = (
                    True if row["label"] == "-1" else False if row["label"] == "1" else None
                )
                yield self._get_peptide_spectrum_match(row)


SageReader = SageTSVReader  # Alias for backwards compatibility


class SageParquetReader(_SageReaderBase):
    """Reader for Sage Parquet results files."""

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with pq.ParquetFile(self.filename) as pq_file:
            for batch in pq_file.iter_batches():
                for row in batch.to_pylist():
                    yield self._get_peptide_spectrum_match(row)


RESCORING_FEATURES = [
    "expmass",
    "calcmass",
    "delta_mass",
    "peptide_len",
    "missed_cleavages",
    "isotope_error",
    "precursor_ppm",
    "fragment_ppm",
    "hyperscore",
    "delta_next",
    "delta_best",
    "delta_rt_model",
    "aligned_rt",
    "predicted_rt",
    "matched_peaks",
    "longest_b",
    "longest_y",
    "longest_y_pct",
    "matched_intensity_pct",
    "scored_candidates",
    "poisson",
    # "ms1_intensity",  # Removed in Sage v0.14
    "ms2_intensity",
]
