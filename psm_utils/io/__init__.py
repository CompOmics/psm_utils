"""
Parsers for proteomics search results from various search engines.

This module provides a unified interface for reading and writing peptide-spectrum match (PSM)
files from various proteomics search engines and analysis tools. It supports automatic file
type detection and conversion between different formats.

The module includes:

- Reader and writer classes for various PSM file formats
- Automatic file type inference from filename patterns
- File conversion utilities
- Progress tracking for long operations
- Type-safe interfaces with comprehensive error handling

Supported file formats include MaxQuant, MS²PIP, Percolator, mzIdentML, pepXML, and many others.
See the documentation for a complete list of supported formats.

Examples
--------
Read a PSM file with automatic format detection:

>>> from psm_utils.io import read_file
>>> psm_list = read_file("results.tsv")

Convert between file formats:

>>> from psm_utils.io import convert
>>> convert("input.msms", "output.mzid")

Write a PSMList to file:

>>> from psm_utils.io import write_file
>>> write_file(psm_list, "output.tsv")

"""

from __future__ import annotations

import re
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Protocol, TypedDict, runtime_checkable

from rich.progress import track

import psm_utils.io.alphadia as alphadia
import psm_utils.io.cbor as cbor
import psm_utils.io.diann as diann
import psm_utils.io.flashlfq as flashlfq
import psm_utils.io.fragpipe as fragpipe
import psm_utils.io.idxml as idxml
import psm_utils.io.ionbot as ionbot
import psm_utils.io.json as json
import psm_utils.io.maxquant as maxquant
import psm_utils.io.msamanda as msamanda
import psm_utils.io.mzid as mzid
import psm_utils.io.parquet as parquet
import psm_utils.io.peptide_record as peptide_record
import psm_utils.io.pepxml as pepxml
import psm_utils.io.percolator as percolator
import psm_utils.io.proteome_discoverer as proteome_discoverer
import psm_utils.io.proteoscape as proteoscape
import psm_utils.io.sage as sage
import psm_utils.io.tsv as tsv
import psm_utils.io.xtandem as xtandem
from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class FileType(TypedDict):
    """Type definition for filetype properties."""

    reader: type[ReaderBase] | None
    writer: type[WriterBase] | None
    extension: str
    filename_pattern: str


FILETYPES: dict[str, FileType] = {
    "flashlfq": {
        "reader": flashlfq.FlashLFQReader,
        "writer": flashlfq.FlashLFQWriter,
        "extension": ".tsv",
        "filename_pattern": r"^.*\.flashlfq\.tsv$",
    },
    "ionbot": {
        "reader": ionbot.IonbotReader,
        "writer": None,
        "extension": "ionbot.first.csv",
        "filename_pattern": r"^ionbot.first.csv$",
    },
    "idxml": {
        "reader": idxml.IdXMLReader,
        "writer": idxml.IdXMLWriter,
        "extension": ".idXML",
        "filename_pattern": r"^.*\.idxml$",
    },
    "msms": {
        "reader": maxquant.MSMSReader,
        "writer": None,
        "extension": "_msms.txt",
        "filename_pattern": r"^.*msms\.txt$",
    },
    "mzid": {
        "reader": mzid.MzidReader,
        "writer": mzid.MzidWriter,
        "extension": ".mzid",
        "filename_pattern": r"^.*\.(?:(?:mzidentml)|(?:mzid))$",
    },
    "peprec": {
        "reader": peptide_record.PeptideRecordReader,
        "writer": peptide_record.PeptideRecordWriter,
        "extension": ".peprec.txt",
        "filename_pattern": r"(^.*\.peprec(?:\.txt)?$)|(?:^peprec\.txt$)",
    },
    "pepxml": {
        "reader": pepxml.PepXMLReader,
        "writer": None,
        "extension": ".pepxml",
        "filename_pattern": r"^.*\.pep\.?xml$",
    },
    "percolator": {
        "reader": percolator.PercolatorTabReader,
        "writer": percolator.PercolatorTabWriter,
        "extension": ".percolator.txt",
        "filename_pattern": r"^.*\.(?:(?:pin)|(?:pout))$",
    },
    "proteome_discoverer": {
        "reader": proteome_discoverer.MSFReader,
        "writer": None,
        "extension": ".msf",
        "filename_pattern": r"^.*\.msf$",
    },
    "proteoscape": {
        "reader": proteoscape.ProteoScapeReader,
        "writer": None,
        "extension": ".parquet",
        "filename_pattern": r"^.*\.candidates\.parquet$",
    },
    "xtandem": {
        "reader": xtandem.XTandemReader,
        "writer": None,
        "extension": ".t.xml",
        "filename_pattern": r"^.*\.t\.xml$",
    },
    "msamanda": {
        "reader": msamanda.MSAmandaReader,
        "writer": None,
        "extension": ".csv",
        "filename_pattern": r"^.*(?:_|\.)msamanda.csv$",
    },
    "sage_tsv": {
        "reader": sage.SageTSVReader,
        "writer": None,
        "extension": ".tsv",
        "filename_pattern": r"^.*(?:_|\.)sage.tsv$",
    },
    "sage_parquet": {
        "reader": sage.SageParquetReader,
        "writer": None,
        "extension": ".parquet",
        "filename_pattern": r"^.*(?:_|\.)sage.parquet$",
    },
    "fragpipe": {
        "reader": fragpipe.FragPipeReader,
        "writer": None,
        "extension": ".tsv",
        "filename_pattern": r"^.*(?:_|\.)?psm\.tsv$",
    },
    "alphadia": {
        "reader": alphadia.AlphaDIAReader,
        "writer": None,
        "extension": ".tsv",
        "filename_pattern": r"^.*(?:_|\.)?precursors\.tsv$",
    },
    "diann": {
        "reader": diann.DIANNTSVReader,
        "writer": None,
        "extension": ".tsv",
        "filename_pattern": r"^.*(?:_|\.)?diann\.tsv$",
    },
    "parquet": {  # List after more specific Parquet patterns to avoid matching conflicts
        "reader": parquet.ParquetReader,
        "writer": parquet.ParquetWriter,
        "extension": ".parquet",
        "filename_pattern": r"^.*\.parquet$",
    },
    "json": {
        "reader": json.JSONReader,
        "writer": json.JSONWriter,
        "extension": ".json",
        "filename_pattern": r"^.*\.json$",
    },
    "cbor": {
        "reader": cbor.CBORReader,
        "writer": cbor.CBORWriter,
        "extension": ".cbor",
        "filename_pattern": r"^.*\.cbor$",
    },
    "tsv": {  # List after more specific TSV patterns to avoid matching conflicts
        "reader": tsv.TSVReader,
        "writer": tsv.TSVWriter,
        "extension": ".tsv",
        "filename_pattern": r"^.*\.tsv$",
    },
}

FILETYPES["sage"] = FILETYPES["sage_tsv"]  # Alias for backwards compatibility

# Type-annotated lookup dictionaries for readers and writers
READERS: dict[str, type[ReaderBase]] = {
    k: v["reader"] for k, v in FILETYPES.items() if v["reader"]
}
WRITERS: dict[str, type[WriterBase]] = {
    k: v["writer"] for k, v in FILETYPES.items() if v["writer"]
}


@runtime_checkable
class _SupportsStr(Protocol):
    """Protocol to check if an object supports string conversion."""

    def __str__(self) -> str: ...


def _infer_filetype(filename: _SupportsStr) -> str:
    """Infer filetype from filename using pattern matching."""
    for filetype, properties in FILETYPES.items():
        if re.fullmatch(properties["filename_pattern"], str(filename), flags=re.IGNORECASE):
            return filetype
    else:
        raise PSMUtilsIOException("Could not infer filetype.")


def _supports_write_psm(writer: type[WriterBase]) -> bool:
    """Check if writer supports write_psm method."""
    with NamedTemporaryFile(delete=False) as temp_file:
        temp_file.close()
        Path(temp_file.name).unlink()
        example_psm = PSM(peptidoform="ACDE", spectrum_id="0")
        try:
            with writer(temp_file.name, example_psm=example_psm) as writer_instance:
                writer_instance.write_psm(example_psm)
        except NotImplementedError:
            supports_write_psm = False
        except AttributeError:  # `None` is not valid PSM
            supports_write_psm = True
        else:
            supports_write_psm = True
            Path(temp_file.name).unlink()
        return supports_write_psm


def read_file(filename: str | Path, *args, filetype: str = "infer", **kwargs) -> PSMList:
    """
    Read PSM file into :py:class:`~psm_utils.psmlist.PSMList`.

    Parameters
    ----------
    filename
        Path to the PSM file to read.
    filetype
        File type specification. Can be any PSM file type with read support or "infer" to
        automatically detect from filename pattern. See documentation for supported file formats.
    *args
        Additional positional arguments passed to the PSM file reader.
    **kwargs
        Additional keyword arguments passed to the PSM file reader.

    Returns
    -------
    List of PSM objects parsed from the input file.

    Raises
    ------
    PSMUtilsIOException
        If filetype cannot be inferred or if the specified filetype is
        unknown or not supported for reading.

    """
    if filetype == "infer":
        filetype = _infer_filetype(filename)
    try:
        reader_cls = READERS[filetype]
    except KeyError as e:
        raise PSMUtilsIOException(
            f"Filetype '{filetype}' unknown or not supported for reading."
        ) from e
    reader = reader_cls(filename, *args, **kwargs)
    psm_list = reader.read_file()
    return psm_list


def write_file(
    psm_list: PSMList,
    filename: str | Path,
    *args,
    filetype: str = "infer",
    show_progressbar: bool = False,
    **kwargs,
) -> None:
    """
    Write :py:class:`~psm_utils.psmlist.PSMList` to PSM file.

    Parameters
    ----------
    psm_list
        List of PSM objects to be written to file.
    filename
        Path to the output file.
    filetype
        File type specification. Can be any PSM file type with write support or "infer" to
        automatically detect from filename pattern. See documentation for supported file formats.
    show_progressbar
        Whether to display a progress bar during the writing process.
    *args
        Additional positional arguments passed to the PSM file writer.
    **kwargs
        Additional keyword arguments passed to the PSM file writer.

    Raises
    ------
    PSMUtilsIOException
        If filetype cannot be inferred or if the specified filetype is
        unknown or not supported for writing.
    IndexError
        If psm_list is empty and cannot provide an example PSM.

    """
    if filetype == "infer":
        filetype = _infer_filetype(filename)
    try:
        writer_cls = WRITERS[filetype]
    except KeyError as e:
        raise PSMUtilsIOException(
            f"Filetype {filetype} unknown or not supported for writing."
        ) from e

    # Remove file if already exists to avoid appending:
    if Path(filename).is_file():
        Path(filename).unlink()

    # Get example PSM, instantiate writer, write
    example_psm = psm_list[0]
    with writer_cls(
        filename,
        *args,
        example_psm=example_psm,
        mode="write",
        show_progressbar=show_progressbar,
        **kwargs,
    ) as writer:
        writer.write_file(psm_list)


def convert(
    input_filename: str | Path,
    output_filename: str | Path,
    input_filetype: str = "infer",
    output_filetype: str = "infer",
    show_progressbar: bool = False,
) -> None:
    """
    Convert a PSM file from one format into another.

    Parameters
    ----------
    input_filename
        Path to the input PSM file.
    output_filename
        Path to the output PSM file.
    input_filetype
        Input file type specification. Can be any PSM file type with read support
        or "infer" to automatically detect from filename pattern.
        See documentation for supported file formats.
    output_filetype
        Output file type specification. Can be any PSM file type with write support
        or "infer" to automatically detect from filename pattern.
        See documentation for supported file formats.
    show_progressbar
        Whether to display a progress bar during the conversion process.

    Raises
    ------
    PSMUtilsIOException
        If input or output filetypes cannot be inferred, if the specified filetypes are
        unknown or not supported, or if the input file is empty.
    KeyError
        If the specified filetype is not found in READERS or WRITERS dictionaries.

    Examples
    --------
    Convert a MaxQuant msms.txt file to a MS²PIP peprec file, while inferring
    the applicable file types from the file extensions:

    >>> from psm_utils.io import convert
    >>> convert("msms.txt", "filename_out.peprec")

    Convert a MaxQuant msms.txt file to a MS²PIP peprec file, while explicitly
    specifying both file types:

    >>> convert(
    ...     "filename_in.msms",
    ...     "filename_out.peprec",
    ...     input_filetype="msms",
    ...     output_filetype="peprec"
    ... )

    Notes
    -----
    Filetypes can only be inferred for select specific file names and/or extensions, such as
    ``msms.txt`` or ``*.peprec``.

    """
    # If needed, infer input and output filetypes
    if input_filetype == "infer":
        input_filetype = _infer_filetype(input_filename)
    if output_filetype == "infer":
        output_filetype = _infer_filetype(output_filename)

    try:
        reader_cls = READERS[input_filetype]
        writer_cls = WRITERS[output_filetype]
    except KeyError as e:
        raise PSMUtilsIOException(f"Filetype '{e.args[0]}' unknown or not supported.") from e

    # Remove file if already exists to avoid appending:
    if Path(output_filename).is_file():
        Path(output_filename).unlink()

    reader = reader_cls(input_filename)

    if _supports_write_psm(writer_cls):
        # Setup iterator, potentially with indeterminate progress bar
        if show_progressbar:
            # Use indeterminate progress tracking for lazy evaluation
            iterator = track(reader, description="[green]Converting file")
        else:
            iterator = reader

        # Get example PSM and instantiate writer
        for psm in reader:
            example_psm = psm
            break
        else:
            raise PSMUtilsIOException("Input file is empty or does not contain valid PSMs.")

        writer = writer_cls(output_filename, example_psm=example_psm, mode="write")

        # Convert
        with writer:
            for psm in iterator:
                writer.write_psm(psm)

    # First read full PSM list, then write file at once
    elif writer_cls == mzid.MzidWriter:
        writer = writer_cls(output_filename, show_progressbar=show_progressbar)
        writer.write_file(reader.read_file())
    else:
        writer = writer_cls(output_filename)
        writer.write_file(reader.read_file())
