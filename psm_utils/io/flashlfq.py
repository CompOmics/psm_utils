"""
Reader and writer for the FlashLFQ generic TSV format.

See the `FlashLFQ documentation <https://github.com/smith-chem-wisc/FlashLFQ/wiki/Identification-Input-Formats>`_
for more information on the format.

Notes
-----
- The FlashLFQ format does not contain the actual spectrum identifier. When reading a FlashLFQ
  file, the spectrum identifier is set to the row number in the file.
- The FlashLFQ format does not contain the precursor m/z, but the theoretical monoisotopic mass.
  This value is not read into the PSM object, but can be calculated from the peptidoform.
- To read from a FlashLFQ file, the ``Full Sequence`` column is expected to contain a ProForma v2
  compatible peptidoform notation.

"""

from __future__ import annotations

import csv
import logging
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()

LOGGER = logging.getLogger(__name__)


class FlashLFQReader(ReaderBase):
    """Reader for FlashLFQ TSV format."""

    required_columns: list[str] = ["Full Sequence", "Precursor Charge"]

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as open_file:
            reader = csv.DictReader(open_file, delimiter="\t")
            if not reader.fieldnames:
                raise PSMUtilsIOException(
                    f"FlashLFQ TSV file '{self.filename}' is empty or has no valid header."
                )
            if not all(col in reader.fieldnames for col in self.required_columns):
                raise PSMUtilsIOException(
                    f"FlashLFQ TSV file must contain the following columns: {self.required_columns}"
                )
            for i, row in enumerate(reader):
                yield self._parse_entry(row, spectrum_id=str(i))

    def _parse_entry(self, entry: dict[str, Any], spectrum_id: str) -> PSM:
        """Parse single FlashLFQ TSV entry to :py:class:`~psm_utils.psm.PSM`."""
        # Replace empty strings with None
        entry = {k: v if v else None for k, v in entry.items()}

        # Parse entry
        return PSM(
            peptidoform=f"{entry['Full Sequence']}/{entry['Precursor Charge']}",
            spectrum_id=spectrum_id,
            run=entry.get("File Name"),
            retention_time=entry.get("Scan Retention Time"),
            protein_list=self._parse_protein_list(entry.get("Protein Accession")),
        )

    @staticmethod
    def _parse_protein_list(protein_accession: str | None) -> list[str]:
        """Parse protein list string to list of protein accessions."""
        if not protein_accession:
            return []
        elif ";" in protein_accession:  # Docs define separator as semicolon
            return protein_accession.split(";")
        elif "|" in protein_accession:  # Example file uses pipe
            return protein_accession.split("|")
        else:
            return [protein_accession]  # Single protein


class FlashLFQWriter(WriterBase):
    """Reader for FlashLFQ TSV format."""

    _default_fieldnames: list[str] = [
        "File Name",
        "Base Sequence",
        "Full Sequence",
        "Peptide Monoisotopic Mass",
        "Scan Retention Time",
        "Precursor Charge",
        "Protein Accession",
    ]

    def __init__(
        self,
        filename: str | Path,
        *args: Any,
        fdr_threshold: float = 0.01,
        only_targets: bool = True,
        **kwargs: Any,
    ) -> None:
        """
        Reader for psm_utils TSV format.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments passed to the base class.
        fdr_threshold
            FDR threshold for filtering PSMs.
        only_targets
            If True, only target PSMs are written to file. If False, both target and decoy PSMs
            are written.
        **kwargs
            Additional keyword arguments passed to the base class.

        """
        super().__init__(filename, *args, **kwargs)

        self.fdr_threshold: float = fdr_threshold
        self.only_targets: bool = only_targets

        self._open_file: Any = None
        self._writer: Any = None
        self.fieldnames: list[str] | None = None

    def __enter__(self) -> FlashLFQWriter:
        """Open file for writing and return self."""
        if Path(self.filename).is_file():
            # Get fieldnames from existing file
            with open(self.filename) as open_file:
                # Get fieldnames
                self.fieldnames = open_file.readline().strip().split("\t")
            mode: str = "at"
        else:
            # Set default fieldnames; avoiding mutation of class variable
            self.fieldnames = self._default_fieldnames[:]
            mode = "wt"

        # Open file and writer
        self._open_file = open(self.filename, mode, newline="")
        self._writer = csv.DictWriter(
            self._open_file,
            fieldnames=self.fieldnames,
            extrasaction="ignore",
            delimiter="\t",
        )

        if mode == "wt":
            self._writer.writeheader()

        return self

    def __exit__(self, *args: Any, **kwargs: Any) -> None:
        """Close file and writer."""
        if self._open_file is not None:
            self._open_file.close()
        self._open_file = None
        self._writer = None

    def write_psm(self, psm: PSM) -> None:
        """
        Write a single PSM to new or existing PSM file.

        Parameters
        ----------
        psm
            PSM object to write.

        """
        if psm.qvalue and psm.qvalue > self.fdr_threshold:
            return
        if self.only_targets and psm.is_decoy:
            return

        entry = self._psm_to_entry(psm)
        try:
            self._writer.writerow(entry)  # type: ignore[union-attr]
        except AttributeError as e:
            raise PSMUtilsIOException(
                f"`write_psm` method can only be called if `{self.__class__.__qualname__}`"
                "is opened in context (i.e., using the `with` statement)."
            ) from e

    def write_file(self, psm_list: PSMList) -> None:
        """
        Write an entire PSMList to a new PSM file.

        Parameters
        ----------
        psm_list
            PSMList object to write to file.

        """
        # Filter out decoys
        if self.only_targets:
            # Accept both None and False
            target_mask = np.array([not psm.is_decoy for psm in psm_list], dtype=bool)
            LOGGER.debug(f"Skipping {(~target_mask).sum()} decoy PSMs for FlashLFQ file.")
        else:
            target_mask = np.ones(len(psm_list), dtype=bool)

        # Filter out PSMs above FDR threshold
        if any(psm.qvalue is None for psm in psm_list):
            LOGGER.warning(
                "Not all PSMs have a q-value. Skipping FDR filtering for FlashLFQ file."
            )
            fdr_mask: np.ndarray[Any, np.dtype[np.bool_]] = np.ones(len(psm_list), dtype=bool)
        else:
            fdr_mask = psm_list["qvalue"] <= self.fdr_threshold
        filtered_by_fdr: int = (~fdr_mask & target_mask).sum()
        LOGGER.debug(f"Skipping {filtered_by_fdr} PSMs above FDR threshold for FlashLFQ file.")

        filtered_psm_list: PSMList = psm_list[target_mask & fdr_mask]

        with open(self.filename, "w", newline="") as f:
            if not self.fieldnames:
                # Set default fieldnames; avoiding mutation of class variable
                self.fieldnames = self._default_fieldnames[:]
            writer = csv.DictWriter(
                f, fieldnames=self.fieldnames, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            for psm in filtered_psm_list:
                writer.writerow(self._psm_to_entry(psm))

    @staticmethod
    def _psm_to_entry(psm: PSM) -> dict[str, Any]:
        """Convert :py:class:`~psm_utils.psm.PSM` to FlashLFQ TSV entry."""
        return {
            "File Name": psm.run,
            "Base Sequence": psm.peptidoform.sequence,
            "Full Sequence": psm.peptidoform.modified_sequence,
            "Peptide Monoisotopic Mass": f"{psm.peptidoform.theoretical_mass:.6f}",
            "Scan Retention Time": psm.retention_time,
            "Precursor Charge": psm.peptidoform.precursor_charge,
            "Protein Accession": ";".join(psm.protein_list) if psm.protein_list else None,
        }
