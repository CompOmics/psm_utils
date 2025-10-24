"""
Interface with ionbot PSM files.

Currently only supports the ionbot.first.csv files.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterator
from pathlib import Path

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM

set_csv_field_size_limit()

REQUIRED_COLUMNS = [
    "database_peptide",
    "modifications",
    "charge",
    "spectrum_title",
    "spectrum_file",
    "proteins",
    "observed_retention_time",
    "database",
    "psm_score",
    "q-value",
    "PEP",
]


class IonbotReader(ReaderBase):
    """Reader for ionbot PSM files."""

    def __init__(
        self,
        filename: str | Path,
        *args,
        **kwargs,
    ) -> None:
        """
        Reader for ``ionbot.first.csv`` PSM files.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments passed to parent class.
        **kwargs
            Additional keyword arguments passed to parent class.

        Examples
        --------
        IonbotReader supports iteration:

        >>> from psm_utils.io.ionbot import IonbotReader
        >>> for psm in IonbotReader("ionbot.first.csv"):
        ...     print(psm.peptidoform.proforma)
        ACDEK
        AC[Carbamidomethyl]DEFGR
        [Acetyl]-AC[Carbamidomethyl]DEFGHIK

        Or a full file can be read at once into a :py:class:`psm_utils.psm_list.PSMList`
        object:

        >>> ionbot_reader = IonbotReader("ionbot.first.csv")
        >>> psm_list = ionbot_reader.read_file()

        """
        super().__init__(filename, *args, **kwargs)

    def __iter__(self) -> Iterator[PSM]:
        """
        Iterate over file and return PSMs one-by-one.

        Yields
        ------
        PSM
            Individual PSM objects from the ionbot CSV file.

        Raises
        ------
        FileNotFoundError
            If the specified file does not exist.
        csv.Error
            If there are issues reading the CSV file.
        InvalidIonbotModificationError
            If modification parsing fails.

        """
        with open(self.filename) as open_file:
            reader = csv.DictReader(open_file, delimiter=",")
            for row in reader:
                yield self._get_peptide_spectrum_match(row)

    def _get_peptide_spectrum_match(self, psm_dict: dict[str, str]) -> PSM:
        """Convert a dictionary row from ionbot CSV to a PSM object."""
        try:
            return PSM(
                peptidoform=self._parse_peptidoform(
                    psm_dict["matched_peptide"],
                    psm_dict["modifications"],
                    psm_dict["charge"],
                ),
                spectrum_id=psm_dict["spectrum_title"],
                run=psm_dict["spectrum_file"],
                is_decoy=(
                    True
                    if psm_dict["database"] == "D"
                    else False
                    if psm_dict["database"] == "T"
                    else None
                ),
                score=float(psm_dict["psm_score"]),
                precursor_mz=float(psm_dict["m/z"]),
                retention_time=float(psm_dict["observed_retention_time"]),
                protein_list=psm_dict["proteins"].split("||"),
                source="ionbot",
                qvalue=float(psm_dict["q-value"]),
                pep=float(psm_dict["PEP"]),
                provenance_data={"ionbot_filename": str(self.filename)},
                metadata={
                    col: str(psm_dict[col])
                    for col in psm_dict.keys()
                    if col not in REQUIRED_COLUMNS
                },
            )
        except KeyError as e:
            raise PSMUtilsIOException(f"Missing required column in ionbot file: {e}") from e
        except ValueError as e:
            raise PSMUtilsIOException(f"Error parsing numeric value in ionbot file: {e}") from e

    @staticmethod
    def _parse_peptidoform(peptide: str, modifications: str, charge: str | int) -> Peptidoform:
        """Parse peptide, modifications, and charge to Peptidoform."""
        # Split peptide into list of amino acids with termini
        peptide_elements: list[str] = [""] + list(peptide) + [""]

        # Add modifications
        pattern: re.Pattern[str] = re.compile(r"^(?P<U>\[\S*?\])?(?P<mod>.*?)(?P<AA>\[\S*?\])?$")

        if modifications:  # Handle empty modifications string
            mod_parts = modifications.split("|")
            if len(mod_parts) % 2 != 0:
                raise InvalidIonbotModificationError(
                    f"Invalid modification string format: '{modifications}'. "
                    "Expected even number of parts (position|label pairs)."
                )

            for position_str, label in zip(mod_parts[::2], mod_parts[1::2]):
                mod_match = pattern.search(label)
                if not mod_match:
                    raise InvalidIonbotModificationError(
                        f"Invalid modification format '{label}' at position {position_str} in "
                        f"'{modifications}'."
                    )

                try:
                    position = int(position_str)
                except ValueError as e:
                    raise InvalidIonbotModificationError(
                        f"Invalid position '{position_str}' in modifications '{modifications}'"
                    ) from e

                if position < 0 or position >= len(peptide_elements):
                    raise InvalidIonbotModificationError(
                        f"Position {position} out of range for peptide '{peptide}' (length {len(peptide_elements) - 2})"
                    )

                if mod_match.group("U"):
                    parsed_label = "U:" + mod_match.group("U")[1:-1]
                else:
                    parsed_label = mod_match.group("mod")
                peptide_elements[position] += f"[{parsed_label}]"

        # Add terminal modifications
        peptide_elements[0] = peptide_elements[0] + "-" if peptide_elements[0] else ""
        peptide_elements[-1] = "-" + peptide_elements[-1] if peptide_elements[-1] else ""
        proforma_seq = "".join(peptide_elements)

        # Add charge state
        proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)


class InvalidIonbotModificationError(PSMUtilsIOException):
    """
    Exception raised when ionbot modification parsing fails.

    This exception is raised when:
    - Modification format is invalid
    - Position values are out of range
    - Modification string structure is malformed
    """

    pass
