"""Interface to MS Amanda CSV result files."""

from __future__ import annotations

import csv
import logging
import re
from collections.abc import Iterator, Sequence
from itertools import compress
from pathlib import Path

import numpy as np

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PSM, Peptidoform

set_csv_field_size_limit()

LOGGER = logging.getLogger(__name__)


# Minimal set of required columns
REQUIRED_COLUMNS = [
    "Title",
    "Sequence",
    "Modifications",
    "Protein Accessions",
    "Amanda Score",
    "m/z",
    "Charge",
    "RT",
    "Filename",
    "Id",
]

# Potential rescoring features
RESCORING_FEATURES = [
    "number of missed cleavages",
    "number of residues",
    "number of considered fragment ions",
    "delta M",
    "avg MS2 error[ppm]",
    "assigned intensity fraction",
    "binom score",
    "Weighted Probability",
    "Nr of matched peaks",
]


class MSAmandaReader(ReaderBase):
    """Reader for MS Amanda CSV result files."""

    def __init__(self, filename: str | Path, *args, **kwargs) -> None:
        """
        Initialize reader for MS Amanda CSV result files.

        Parameters
        ----------
        filename
            Path to the MS Amanda CSV file.
        *args
            Additional positional arguments passed to parent class.
        **kwargs
            Additional keyword arguments passed to parent class.

        """
        super().__init__(filename, *args, **kwargs)
        self._present_columns: list[str] = REQUIRED_COLUMNS.copy()
        self._rescoring_feature_columns: list[str] = []
        self._metadata_columns: list[str] = []
        self._has_rank_column: bool | None = None

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as open_file:
            if not next(open_file).startswith("#"):
                open_file.seek(0)
            reader = csv.DictReader(open_file, delimiter="\t")
            self._evaluate_columns(reader.fieldnames)
            for psm_dict in reader:
                yield self._get_peptide_spectrum_match(psm_dict)

    def _evaluate_columns(self, columns: Sequence[str] | None) -> None:
        """Evaluate and validate columns from MS Amanda file header."""
        if columns is None:
            raise MSAmandaParsingError("MS Amanda file does not contain any columns.")

        # Check if required columns are present
        column_check = [True if col in columns else False for col in REQUIRED_COLUMNS]
        if not all(column_check):
            missing = list(compress(REQUIRED_COLUMNS, list(~np.array(column_check))))
            raise MSAmandaParsingError(f"Missing columns: {missing}")

        # Check if rank is present
        if "Rank" in columns:
            self._has_rank_column = True
            self._present_columns.append("Rank")

        # Get list of present rescoring features
        self._rescoring_feature_columns = [col for col in RESCORING_FEATURES if col in columns]

        # Add remaining columns to metadata
        self._metadata_columns = [
            col
            for col in columns
            if col not in self._present_columns + self._rescoring_feature_columns
        ]

    def _get_peptide_spectrum_match(self, psm_dict: dict[str, str]) -> PSM:
        """Return a PSM object from MS Amanda CSV PSM file."""
        psm = PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["Sequence"], psm_dict["Modifications"], psm_dict["Charge"]
            ),
            spectrum_id=psm_dict["Title"],
            run=Path(psm_dict["Filename"]).with_suffix("").name,
            is_decoy=psm_dict["Protein Accessions"].startswith("REV_"),
            score=float(psm_dict["Amanda Score"]),
            precursor_mz=float(psm_dict["m/z"]),
            retention_time=float(psm_dict["RT"]),
            protein_list=psm_dict["Protein Accessions"].split(";"),
            source="MSAmanda",
            provenance_data=({"MSAmanda_filename": str(self.filename)}),
            rescoring_features={
                col: str(value)
                for col, value in psm_dict.items()
                if col in self._rescoring_feature_columns
            },
            metadata={
                col: str(value) for col, value in psm_dict.items() if col in self._metadata_columns
            },
        )
        if self._has_rank_column:
            psm["rank"] = psm_dict["Rank"]
        return psm

    @staticmethod
    def _parse_peptidoform(seq: str, modifications: str, charge: str) -> Peptidoform:
        """Parse MSAmanda sequence, modifications and charge to ProForma sequence."""
        peptide = [""] + [aa.upper() for aa in seq] + [""]
        pattern = re.compile(
            r"(?:(?:(?P<site>[A-Z])(?P<loc>\d+))|(?P<term>[CN]-Term))\((?P<mod_name>[^|()]+)\|(?P<mz>[-0-9.]+)\|(?P<type>variable|fixed)\);?"
        )

        for match in pattern.finditer(modifications):
            if match.group("term") == "N-Term":
                peptide[0] = peptide[0] + f"[{match.group('mod_name')}]"
            elif match.group("term") == "C-Term":
                peptide[-1] = peptide[-1] + f"[{match.group('mod_name')}]"
            if match.group("loc") is not None:
                peptide[int(match.group("loc"))] = (
                    peptide[int(match.group("loc"))] + f"[{match.group('mod_name')}]"
                )

        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""
        proforma_seq = "".join(peptide)

        return Peptidoform(proforma_seq + f"/{charge}")


class MSAmandaParsingError(PSMUtilsIOException):
    """Error in parsing MS Amanda CSV file."""

    pass
