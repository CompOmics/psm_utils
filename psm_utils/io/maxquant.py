"""Interface to MaxQuant msms.txt PSM files."""

from __future__ import annotations

import csv
import logging
import re
from collections.abc import Iterator, Sequence
from itertools import compress
from pathlib import Path

import numpy as np

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM

set_csv_field_size_limit()

logger = logging.getLogger(__name__)

MSMS_REQUIRED_COLUMNS = [
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


class MSMSReader(ReaderBase):
    """Reader for MaxQuant msms.txt PSM files."""

    def __init__(
        self,
        filename: str | Path,
        *args,
        **kwargs,
    ) -> None:
        """
        Initialize reader for MaxQuant msms.txt PSM files.

        Parameters
        ----------
        filename
            Path to the MaxQuant msms.txt PSM file.
        *args
            Additional positional arguments passed to parent class.
        **kwargs
            Additional keyword arguments passed to parent class.

        Examples
        --------
        MSMSReader supports iteration:

        >>> from psm_utils.io.maxquant import MSMSReader
        >>> for psm in MSMSReader("msms.txt"):
        ...     print(psm.peptidoform.proforma)
        WFEELSK
        NDVPLVGGK
        GANLGEMTNAGIPVPPGFC[+57.022]VTAEAYK
        ...

        Or a full file can be read at once into a :py:class:`~psm_utils.psm_list.PSMList`
        object:

        >>> reader = MSMSReader("msms.txt")
        >>> psm_list = reader.read_file()

        """
        super().__init__(filename, *args, **kwargs)
        self._validate_msms()

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")
            for psm_dict in reader:
                psm = self._get_peptide_spectrum_match(psm_dict)
                yield psm

    def _validate_msms(self) -> None:
        """Validate that the msms.txt file contains required columns."""
        with open(self.filename) as msms_file:
            msms_reader = csv.DictReader(msms_file, delimiter="\t")
            self._evaluate_columns(msms_reader.fieldnames)

    @staticmethod
    def _evaluate_columns(columns: Sequence[str] | None) -> None:
        """Case insensitive column evaluation for msms file."""
        if columns is None:
            raise MSMSParsingError("MSMS file does not contain any columns.")
        columns = list(map(lambda col: col.lower(), columns))
        column_check = [True if col.lower() in columns else False for col in MSMS_REQUIRED_COLUMNS]
        if not all(column_check):
            raise MSMSParsingError(
                f"Missing columns: {list(compress(MSMS_REQUIRED_COLUMNS, list(~np.array(column_check))))}"
            )

    def _get_peptide_spectrum_match(self, psm_dict: dict[str, str]) -> PSM:
        """Return a PSM object from MaxQuant msms.txt PSM file."""
        psm = PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["Modified sequence"], int(psm_dict["Charge"])
            ),
            spectrum_id=psm_dict["Scan number"],
            run=psm_dict["Raw file"],
            is_decoy=psm_dict["Reverse"] == "+",
            score=float(psm_dict["Score"]),
            precursor_mz=float(psm_dict["m/z"]),
            retention_time=float(psm_dict["Retention time"]),
            protein_list=psm_dict["Proteins"].split(";") if psm_dict["Proteins"] else None,
            pep=psm_dict["PEP"],
            source="msms",
            provenance_data=({"msms_filename": str(self.filename)}),
            metadata={
                col: str(psm_dict[col])
                for col in psm_dict.keys()
                if col not in MSMS_REQUIRED_COLUMNS
            },
        )
        return psm

    @staticmethod
    def _parse_peptidoform(modified_seq: str, charge: int) -> Peptidoform:
        """Parse modified sequence to Peptidoform."""
        # pattern to match open and closed round brackets
        pattern = re.compile(r"\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)")
        modified_seq = modified_seq.strip("_")
        modified_seq_len = len(modified_seq)

        for match in pattern.finditer(modified_seq):
            # take match object string, remove leading and trailing bracket
            # and escape remaining brackets
            se_mod_string = match[1].replace("(", r"\(").replace(")", r"\)")

            # if N-term mod
            if match.start() == 0:
                modified_seq = re.sub(rf"\({se_mod_string}\)", rf"[{match[1]}]-", modified_seq)

            # if C-term mod
            elif match.end() == modified_seq_len:
                modified_seq = re.sub(rf"\({se_mod_string}\)", rf"-[{match[1]}]", modified_seq)

            # if modification on amino acid
            else:
                modified_seq = re.sub(rf"\({se_mod_string}\)", rf"[{match[1]}]", modified_seq)

        modified_seq += f"/{charge}"

        return Peptidoform(modified_seq)


class MSMSParsingError(PSMUtilsException):
    """Error while parsing MaxQuant msms.txt PSM file."""
