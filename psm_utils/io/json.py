"""
Reader and writer for a simple, lossless psm_utils JSON format.

Similar to the :py:mod:`psm_utils.io.tsv` and :py:mod:`psm_utils.io.parquet` modules,
this module provides a reader and writer for :py:class:`~psm_utils.psm_list.PSMList`
objects in a lossless manner using JSON format. JSON provides human-readable output
and is widely compatible across platforms and programming languages.

The JSON format stores PSMs as an array of objects, where each object represents a PSM
with its attributes. Peptidoforms are written in the `HUPO-PSI ProForma 2.0
<https://psidev.info/proforma>`_ notation. Fields that are not set (i.e., have a value of None)
are omitted from the JSON output to reduce file size.

"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

from pydantic import ValidationError

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)


class JSONReader(ReaderBase):
    """Reader for psm_utils JSON format."""

    def __init__(self, filename: str | Path, *args, **kwargs):
        """
        Reader for psm_utils JSON format.

        Parameters
        ----------
        filename: str, Pathlib.Path
            Path to PSM file.

        """
        super().__init__(filename, *args, **kwargs)

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename, "rt", encoding="utf-8") as open_file:
            data = json.load(open_file)

            if not isinstance(data, list):
                raise PSMUtilsIOException("JSON file must contain an array of PSM objects.")

            failed_rows = 0
            for row in data:
                try:
                    yield PSM(**self._parse_entry(row))
                except ValidationError as e:
                    failed_rows += 1
                    logger.warning(f"Could not parse PSM from entry: `{row}`")
                    if failed_rows >= 3:
                        raise PSMUtilsIOException(
                            "Could not parse PSM from three consecutive entries. Verify that the "
                            "file is formatted correctly as a psm_utils JSON file or that the "
                            "correct file type reader is used."
                        ) from e
                else:
                    failed_rows = 0

    @staticmethod
    def _parse_entry(entry: dict) -> dict:
        """Parse single JSON entry to :py:class:`~psm_utils.psm.PSM`."""
        # Ensure dict properties have default values if missing
        if "provenance_data" not in entry:
            entry["provenance_data"] = {}
        if "metadata" not in entry:
            entry["metadata"] = {}
        if "rescoring_features" not in entry:
            entry["rescoring_features"] = {}

        return entry


class JSONWriter(WriterBase):
    """Writer for psm_utils JSON format."""

    def __init__(
        self,
        filename: str | Path,
        indent: Optional[int] = 2,
        *args,
        **kwargs,
    ):
        """
        Writer for psm_utils JSON format.

        Parameters
        ----------
        filename: str, Pathlib.Path
            Path to PSM file.
        indent: int, optional
            Number of spaces for JSON indentation. Set to None for compact output.
            Default is 2.

        """
        super().__init__(filename, *args, **kwargs)
        self.indent = indent
        self._psm_cache = []

    def __enter__(self) -> JSONWriter:
        return self

    def __exit__(self, *args, **kwargs) -> None:
        self._flush()

    def write_psm(self, psm: PSM):
        """Write a single PSM to the JSON file."""
        self._psm_cache.append(self._psm_to_entry(psm))

    def write_file(self, psm_list: PSMList):
        """Write an entire PSMList to the JSON file."""
        for psm in psm_list.psm_list:
            self.write_psm(psm)
        self._flush()

    @staticmethod
    def _psm_to_entry(psm: PSM) -> dict:
        """Convert PSM to a dictionary suitable for JSON serialization."""
        psm_dict = dict(psm)
        # Convert peptidoform to string
        psm_dict["peptidoform"] = str(psm.peptidoform)
        # Remove None values to reduce file size
        psm_dict = {k: v for k, v in psm_dict.items() if v is not None}
        return psm_dict

    def _flush(self):
        """Write the cached PSMs to the JSON file."""
        if not self._psm_cache:
            self._psm_cache = []

        with open(self.filename, "wt", encoding="utf-8") as open_file:
            json.dump(self._psm_cache, open_file, indent=self.indent)

        self._psm_cache = []
