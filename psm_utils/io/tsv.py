"""
Reader and writer for a simple, lossless psm_utils TSV format.

Most PSM file formats will introduce a loss of some information when reading,
writing, or converting with :py:mod:`psm_utils.io` due to differences between file
formats. In contrast, :py:class:`~psm_utils.io.psm_list.PSMList` objects can be written
to — or read from — this simple TSV format without any information loss (with exception
of the free-form :py:attr:`spectrum` attribute).

The format follows basic TSV rules, using tab as delimiter, and supports quoting when
a field contains the delimiter. Peptidoforms are written in the `HUPO-PSI ProForma 2.0
<https://psidev.info/proforma>`_ notation.

Required and optional columns equate to the required and optional attributes of
:py:class:`~psm_utils.psm.PSM`. Dictionary items in
:py:attr:`provenance_data`, :py:attr:`metadata`, and :py:attr:`rescoring_features`
are flattened to separate columns, each with their column names prefixed with
``provenance:``, ``meta:``, and ``rescoring:``, respectively.


**Examples**

.. code-block::
    :caption: Minimal :py:mod:`psm_utils` TSV file

    peptidoform	spectrum_id
    RNVIDKVAK/2	1
    KHLEQHPK/2	2
    ...

.. code-block::
    :caption: Recommended :py:mod:`psm_utils` TSV file, compatible with `HUPO-PSI Universal Spectrum Identifier <https://www.psidev.info/usi>`_

    peptidoform	spectrum_id	run	collection
    VLHPLEGAVVIIFK/2	17555	Adult_Frontalcortex_bRP_Elite_85_f09	PXD000561
    ...

.. code-block::
    :caption: Full :py:mod:`psm_utils` TSV file, converted from a Percolator Tab file

    peptidoform	spectrum_id	run	collection	spectrum	is_decoy	score	precursor_mz	retention_time	protein_list	source	provenance:filename	rescoring:ExpMass	rescoring:CalcMass	rescoring:hyperscore	rescoring:deltaScore	rescoring:frac_ion_b	rescoring:frac_ion_y	rescoring:Mass	rescoring:dM	rescoring:absdM	rescoring:PepLen	rescoring:Charge2	rescoring:Charge3	rescoring:Charge4	rescoring:enzN	rescoring:enzC	rescoring:enzInt
    RNVIDKVAK/2	_3_2_1				False	20.3	1042.64		['DECOY_sp|Q8U0H4_REVERSED|RTCB_PYRFU-tRNA-splicing-ligase-RtcB-OS=Pyrococcus-furiosus...']	percolator	pyro.t.xml.pin	1042.64	1042.64	20.3	6.6	0.444444	0.333333	1042.64	0.0003	0.0003	9	1	0	0	1	0	1
    KHLEQHPK/2	_4_2_1				False	26.5	1016.56		['sp|Q8TZD9|RS15_PYRFU-30S-ribosomal-protein-S15-OS=Pyrococcus-furiosus-(strain-ATCC...']	percolator	pyro.t.xml.pin	1016.56	1016.56	26.5	18.5	0.375	0.75	1016.56	0.001	0.001	8	1	0	0	1	0	0
    ...


"""

from __future__ import annotations

import ast
import csv
import logging
from collections.abc import Iterator
from pathlib import Path
from typing import Any, TextIO

from pydantic import ValidationError

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()

logger = logging.getLogger(__name__)


class TSVReader(ReaderBase):
    """Reader for psm_utils TSV format."""

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as open_file:
            reader = csv.DictReader(open_file, delimiter="\t")
            failed_rows = 0
            for row in reader:
                try:
                    yield PSM(**self._parse_entry(row))
                except ValidationError as e:
                    failed_rows += 1
                    logger.warning(f"Could not parse PSM from row: `{row}`")
                    if failed_rows >= 3:
                        raise PSMUtilsIOException(
                            "Could not parse PSM from three consecutive rows. Verify that the "
                            "file is formatted correctly as a psm_utils TSV file or that the "
                            "correct file type reader is used."
                        ) from e
                else:
                    failed_rows = 0

    @staticmethod
    def _parse_entry(entry: dict[str, str | None]) -> dict[str, Any]:
        """Parse single TSV entry to PSM dict."""
        # Replace empty strings with None
        entry = {k: v if v else None for k, v in entry.items()}

        # Parse protein list
        if "protein_list" in entry and entry["protein_list"]:
            try:
                entry["protein_list"] = ast.literal_eval(entry["protein_list"])
            except ValueError as e:
                raise PSMUtilsIOException(
                    f"Could not parse protein list: `{entry['protein_list']}`"
                ) from e

        # Extract dict properties
        parsed_entry: dict[str, Any] = {}
        provenance_data: dict[str, str | None] = {}
        metadata: dict[str, str | None] = {}
        rescoring_features: dict[str, str | None] = {}
        for k, v in entry.items():
            if k.startswith("provenance:"):
                provenance_data[k[11:]] = v
            elif k.startswith("meta:"):
                metadata[k[5:]] = v
            elif k.startswith("rescoring:"):
                rescoring_features[k[10:]] = v
            else:
                parsed_entry[k] = v

        parsed_entry.update(
            {
                "provenance_data": provenance_data,
                "metadata": metadata,
                "rescoring_features": rescoring_features,
            }
        )

        return parsed_entry


class TSVWriter(WriterBase):
    """Writer for psm_utils TSV format."""

    def __init__(
        self,
        filename: str | Path,
        *args: Any,
        example_psm: PSM | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Writer for psm_utils TSV format.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments passed to the base class.
        example_psm
            Example PSM, required to extract the column names when writing to a new
            file. Should contain all fields that are to be written to the PSM file,
            i.e., all items in the :py:attr:`provenance_data`, :py:attr:`metadata`, and
            :py:attr:`rescoring_features` attributes. In other words, items that are
            not present in the example PSM will not be written to the file, even though
            they are present in other PSMs passed to :py:meth:`write_psm` or
            :py:meth:`write_file`.
        **kwargs
            Additional keyword arguments passed to the base class.

        """
        super().__init__(filename, *args, **kwargs)

        self._open_file: TextIO | None = None
        self._writer: csv.DictWriter[str] | None = None

        if example_psm:
            self.fieldnames: list[str] | None = list(self._psm_to_entry(example_psm).keys())
        else:
            self.fieldnames = None

    def __enter__(self) -> TSVWriter:
        """Enter context manager for file writing."""
        if Path(self.filename).is_file():
            with open(self.filename) as open_file:
                # Get fieldnames
                self.fieldnames = open_file.readline().strip().split("\t")
            self._open_file = open(self.filename, "a", newline="")
            self._writer = csv.DictWriter(
                self._open_file,
                fieldnames=self.fieldnames,
                extrasaction="ignore",
                delimiter="\t",
            )
        else:
            if not self.fieldnames:
                raise ValueError("`example_psm` required when writing to new file.")
            self._open_file = open(self.filename, "w", newline="")
            self._writer = csv.DictWriter(
                self._open_file,
                fieldnames=self.fieldnames,
                extrasaction="ignore",
                delimiter="\t",
            )
            self._writer.writeheader()
        return self

    def __exit__(self, *args: Any, **kwargs: Any) -> None:
        """Exit context manager and clean up file resources."""
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
        if self._writer is None:
            raise PSMUtilsIOException(
                f"`write_psm` method can only be called if `{self.__class__.__qualname__}`"
                "is opened in context (i.e., using the `with` statement)."
            )
        entry = self._psm_to_entry(psm)
        self._writer.writerow(entry)

    def write_file(self, psm_list: PSMList) -> None:
        """
        Write an entire PSMList to a new PSM file.

        Parameters
        ----------
        psm_list
            PSMList object to write to file.

        """
        if not self.fieldnames:
            raise ValueError("`example_psm` required when writing to new file.")
        with open(self.filename, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=self.fieldnames, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            for psm in psm_list:
                writer.writerow(self._psm_to_entry(psm))

    @staticmethod
    def _psm_to_entry(psm: PSM) -> dict[str, Any]:
        """Convert PSM object to dictionary entry for TSV writing."""
        entry = psm.__dict__.copy()

        # Convert Peptidoform to proforma sequence
        entry["peptidoform"] = psm.peptidoform.proforma

        # Drop spectrum
        del entry["spectrum"]

        # Flatten dictionary items
        if entry["provenance_data"]:
            entry.update({"provenance:" + k: v for k, v in entry["provenance_data"].items()})
        if entry["metadata"]:
            entry.update({"meta:" + k: v for k, v in entry["metadata"].items()})
        if entry["rescoring_features"]:
            entry.update({"rescoring:" + k: v for k, v in entry["rescoring_features"].items()})
        del entry["provenance_data"]
        del entry["metadata"]
        del entry["rescoring_features"]

        return entry
