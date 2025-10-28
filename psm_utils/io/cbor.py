"""
Reader and writer for a simple, lossless psm_utils CBOR format.

Similar to the :py:mod:`psm_utils.io.json` module, this module provides a reader and
writer for :py:class:`~psm_utils.psm_list.PSMList` objects in a lossless manner using
CBOR (Concise Binary Object Representation) format. CBOR provides better performance
and smaller file sizes compared to JSON while maintaining similar data structures.

The CBOR format stores PSMs as an array of objects, where each object represents a PSM
with its attributes. Peptidoforms are written in the `HUPO-PSI ProForma 2.0
<https://psidev.info/proforma>`_ notation. Fields that are not set (i.e., have a value of None)
are omitted from the CBOR output to reduce file size.

Note: This module requires the ``cbor2`` package to be installed.

"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from pydantic import ValidationError

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

try:
    import cbor2  # type: ignore[import]

    _has_cbor2 = True
except ImportError:
    _has_cbor2 = False

logger = logging.getLogger(__name__)


class CBORReader(ReaderBase):
    """Reader for psm_utils CBOR format."""

    def __init__(self, filename: str | Path, *args, **kwargs):
        """
        Reader for psm_utils CBOR format.

        Parameters
        ----------
        filename: str, Pathlib.Path
            Path to PSM file.
        *args
            Additional positional arguments passed to the base class.
        **kwargs
            Additional keyword arguments passed to the base class.

        """
        super().__init__(filename, *args, **kwargs)
        if not _has_cbor2:
            raise ImportError(
                "The cbor2 package is required to use the CBOR reader/writer. "
                "Install it with: pip install cbor2"
            )

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename, "rb") as open_file:
            try:
                data = cbor2.load(open_file)  # type: ignore[attr-defined]
            except cbor2.CBORDecodeError as e:  # type: ignore[attr-defined]
                raise PSMUtilsIOException(f"Could not decode CBOR file: {e}") from e

            if not isinstance(data, list):
                raise PSMUtilsIOException("CBOR file must contain an array of PSM objects.")

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
                            "file is formatted correctly as a psm_utils CBOR file or that the "
                            "correct file type reader is used."
                        ) from e
                else:
                    failed_rows = 0

    @staticmethod
    def _parse_entry(entry: dict) -> dict:
        """Parse single CBOR entry to :py:class:`~psm_utils.psm.PSM`."""
        # Ensure dict properties have default values if missing
        if "provenance_data" not in entry:
            entry["provenance_data"] = {}
        if "metadata" not in entry:
            entry["metadata"] = {}
        if "rescoring_features" not in entry:
            entry["rescoring_features"] = {}

        return entry


class CBORWriter(WriterBase):
    """Writer for psm_utils CBOR format."""

    def __init__(
        self,
        filename: str | Path,
        *args,
        **kwargs,
    ):
        """
        Writer for psm_utils CBOR format.

        Parameters
        ----------
        filename: str, Pathlib.Path
            Path to PSM file.
        *args
            Additional positional arguments passed to the base class.
        **kwargs
            Additional keyword arguments passed to the base class.

        """
        super().__init__(filename, *args, **kwargs)
        if not _has_cbor2:
            raise ImportError(
                "The cbor2 package is required to use the CBOR reader/writer. "
                "Install it with: pip install cbor2"
            )
        self._psm_cache: list[dict[str, Any]] = []

    def __enter__(self) -> CBORWriter:
        """Enter context manager."""
        return self

    def __exit__(self, *args, **kwargs) -> None:
        """Exit context manager and flush remaining PSMs to file."""
        self._flush()

    def write_psm(self, psm: PSM):
        """Write a single PSM to the CBOR file."""
        self._psm_cache.append(self._psm_to_entry(psm))

    def write_file(self, psm_list: PSMList):
        """Write an entire PSMList to the CBOR file."""
        for psm in psm_list.psm_list:
            self.write_psm(psm)
        self._flush()

    @staticmethod
    def _psm_to_entry(psm: PSM) -> dict:
        """Convert PSM to a dictionary suitable for CBOR serialization."""
        psm_dict = dict(psm)
        # Convert peptidoform to string
        psm_dict["peptidoform"] = str(psm.peptidoform)
        # Remove None values to reduce file size
        psm_dict = {k: v for k, v in psm_dict.items() if v is not None}
        return psm_dict

    def _flush(self):
        """Write the cached PSMs to the CBOR file."""
        with open(self.filename, "wb") as open_file:
            cbor2.dump(self._psm_cache, open_file)  # type: ignore[attr-defined]

        self._psm_cache = []
