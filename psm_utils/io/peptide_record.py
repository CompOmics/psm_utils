"""
Interface with Peptide Record PSM files.

Peptide Record (or PEPREC) is a legacy PSM file type developed at CompOmics as input
format for `MS²PIP <https://github.com/compomics/ms2pip_c/>`_. It is a simple and
flexible delimited text file where each row represents a single PSM. Required columns
are:

- ``spec_id``: Spectrum identifier; usually the identifier used in the spectrum file.
- ``peptide``: Simple, stripped peptide sequence (e.g., ``ACDE``).
- ``modifications``: Amino acid modifications in a custom format (see below).

Depending on the use case, more columns can be required or optional:

- ``charge``: Peptide precursor charge.
- ``observed_retention_time``: Observed retention time.
- ``predicted_retention_time``: Predicted retention time.
- ``label``: Target/decoy: ``1`` for target PSMs, ``-1`` for decoy PSMs.
- ``score``: Primary search engine score (e.g., the score used for q-value calculation).

Peptide modifications are denoted as a pipe-separated list of pipe-separated
*location → label* pairs for each modification. The location is an integer counted
starting at 1 for the first amino acid. ``0`` is reserved for N-terminal modifications
and ``-1`` for C-terminal modifications. Unmodified peptides can be marked with a hyphen
(``-``). For example:

==================================== ==============================================================================
 PEPREC modification(s)               Explanation
==================================== ==============================================================================
 ``-``                               Unmodified
 ``1|Oxidation``                     ``Oxidation`` on the first amino acid
 ``1|Oxidation|5|Carbamidomethyl``   ``Oxidation`` on the first amino acid and ``Carbamidomethyl`` on the fifth
 ``0|Acetylation``                   ``Acetylation`` on the N-terminus
 ``-1|Amidation``                    ``Amidation`` on the C-terminus
==================================== ==============================================================================


Full PEPREC example:

.. code-block::

    spec_id,modifications,peptide,charge
    peptide1,-,ACDEK,2
    peptide2,2|Carbamidomethyl,ACDEFGR,3
    peptide3,0|Acetyl|2|Carbamidomethyl,ACDEFGHIK,2

.. Attention::
    Labile, unlocalized, and fixed modifications are not encoded in the Peptide Record
    notation. To encode fixed modifications, use
    :py:meth:`~psm_utils.psm_list.PSMList.apply_fixed_modifications` before writing to
    Peptide Record.

"""

from __future__ import annotations

import csv
from collections.abc import Iterator
from pathlib import Path
from typing import Any, TextIO

import pandas as pd
from pydantic import BaseModel, ConfigDict

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()


_REQUIRED_COLUMNS = ["spec_id", "peptide", "modifications"]
_OPTIONAL_COLUMNS = [
    "charge",
    "observed_retention_time",
    "predicted_retention_time",
    "label",
    "score",
]


def _analyze_peprec_file(filename: str | Path) -> tuple[str, list[str]]:
    """Analyze Peptide Record file to determine separator and validate columns."""
    separator = ""
    header: list[str] = []

    with open(filename) as f:
        line = f.readline().strip()
        for sep in ["\t", ",", ";", " "]:
            cols = line.split(sep)
            if all(rc in cols for rc in _REQUIRED_COLUMNS):
                separator = sep
                header = cols
                break
        else:
            raise InvalidPeprecError(
                "Could not infer separator. Please validate the Peptide Record "
                "header and/or the required columns."
            )

    # Validate required columns
    for rc in _REQUIRED_COLUMNS:
        if rc not in header:
            raise InvalidPeprecError(f"Required column missing: `{rc}`")

    return separator, header


class _PeprecEntry(BaseModel):
    """Typed entry for Peptide Record data."""

    spec_id: str
    peptide: str
    modifications: str
    charge: str | None = None
    observed_retention_time: str | None = None
    predicted_retention_time: str | None = None
    label: str | None = None
    score: str | None = None

    model_config = ConfigDict(extra="ignore")


class PeptideRecordReader(ReaderBase):
    """Reader for Peptide Record PSM files."""

    def __init__(
        self,
        filename: str | Path,
        *args,
        **kwargs,
    ) -> None:
        """
        Reader for Peptide Record PSM files.

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
        PeptideRecordReader supports iteration:

        >>> from psm_utils.io.peptide_record import PeptideRecordReader
        >>> for psm in PeptideRecordReader("peprec.txt"):
        ...     print(psm.peptidoform.proforma)
        ACDEK
        AC[Carbamidomethyl]DEFGR
        [Acetyl]-AC[Carbamidomethyl]DEFGHIK

        Or a full file can be read at once into a :py:class:`~psm_utils.psm_list.PSMList` object:

        >>> peprec_reader = PeptideRecordReader("peprec.txt")
        >>> psm_list = peprec_reader.read_file()

        """
        super().__init__(filename, *args, **kwargs)
        self._separator, self._header = _analyze_peprec_file(self.filename)

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as open_file:
            reader = csv.DictReader(open_file, delimiter=self._separator)
            for row in reader:
                entry = _PeprecEntry(**row)
                psm = self._entry_to_psm(entry, filename=self.filename)
                yield psm

    @staticmethod
    def _entry_to_psm(entry: _PeprecEntry, filename: str | Path) -> PSM:
        """Parse single Peptide Record entry to PSM."""
        # Parse sequence and modifications
        charge = int(entry.charge) if entry.charge else None
        proforma = peprec_to_proforma(entry.peptide, entry.modifications, charge)

        # Parse decoy label
        if entry.label:
            is_decoy_map = {"-1": True, "1": False}
            try:
                is_decoy = is_decoy_map[entry.label]
            except (ValueError, KeyError) as e:
                raise InvalidPeprecError(
                    f"Could not parse value for `label` {entry.label}. Should be `1` or `-1`."
                ) from e
        else:
            is_decoy = None

        return PSM(
            peptidoform=proforma,
            spectrum_id=entry.spec_id,
            is_decoy=is_decoy,
            retention_time=entry.observed_retention_time,
            score=entry.score,
            source="PeptideRecord",
            provenance_data={"peprec_filename": str(filename)},
        )


class PeptideRecordWriter(WriterBase):
    """Writer for Peptide Record PSM files."""

    def __init__(self, filename: str | Path, *args, **kwargs) -> None:
        """
        Writer for Peptide Record PSM files.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments passed to parent class.
        **kwargs
            Additional keyword arguments passed to parent class.

        """
        super().__init__(filename, *args, **kwargs)
        self._open_file: TextIO | None = None
        self._writer: csv.DictWriter | None = None

    def __enter__(self) -> PeptideRecordWriter:
        """Open file for writing and prepare CSV writer."""
        # If file exists, analyze it to determine separator and header
        if Path(self.filename).is_file():
            separator, header = _analyze_peprec_file(self.filename)
            self._open_file = open(self.filename, "a", newline="")
            self._writer = csv.DictWriter(
                self._open_file,
                fieldnames=header,
                extrasaction="ignore",
                delimiter=separator,
            )
        else:
            self._open_file = open(self.filename, "w", newline="")
            self._writer = csv.DictWriter(
                self._open_file,
                fieldnames=_REQUIRED_COLUMNS + _OPTIONAL_COLUMNS,
                extrasaction="ignore",
                delimiter=" ",
            )
            self._writer.writeheader()
        return self

    def __exit__(self, *args, **kwargs) -> None:
        """Close file when exiting context."""
        if self._open_file is not None:
            self._open_file.close()
        self._open_file = None
        self._writer = None

    @staticmethod
    def _psm_to_entry(psm: PSM) -> dict[str, Any]:
        """Convert PSM to Peptide Record entry dictionary."""
        sequence, modifications, charge = proforma_to_peprec(psm.peptidoform)
        return {
            "spec_id": psm.spectrum_id,
            "peptide": sequence,
            "modifications": modifications,
            "charge": charge,
            "label": None if psm.is_decoy is None else -1 if psm.is_decoy else 1,
            "observed_retention_time": psm.retention_time,
            "score": psm.score,
        }

    def write_psm(self, psm: PSM) -> None:
        """
        Write a single PSM to new or existing Peptide Record PSM file.

        Parameters
        ----------
        psm
            PSM object to write.

        Examples
        --------
        To write single PSMs to a file, :py:class:`PeptideRecordWriter` must be opened as
        a context manager. Then, within the context, :py:func:`write_psm` can be called:

        >>> with PeptideRecordWriter("peprec.txt") as writer:
        >>>     writer.write_psm(psm)

        """
        if self._writer is None:
            raise PSMUtilsIOException(
                f"`write_psm` method can only be called if `{self.__class__.__qualname__}` "
                "is opened in context (i.e., using the `with` statement)."
            )
        entry = self._psm_to_entry(psm)
        try:
            self._writer.writerow(entry)
        except AttributeError as e:
            raise PSMUtilsIOException(
                f"`write_psm` method can only be called if `{self.__class__.__qualname__}` "
                "is opened in context (i.e., using the `with` statement)."
            ) from e

    def write_file(self, psm_list: PSMList) -> None:
        """
        Write an entire PSMList to a new Peptide Record PSM file.

        Parameters
        ----------
        psm_list
            PSMList object to write to file.

        Examples
        --------
        >>> writer = PeptideRecordWriter("peprec.txt")
        >>> writer.write_file(psm_list)

        """
        with open(self.filename, "w", newline="") as f:
            fieldnames = _REQUIRED_COLUMNS + _OPTIONAL_COLUMNS
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=" ")
            writer.writeheader()
            for psm in psm_list:
                writer.writerow(self._psm_to_entry(psm))


def peprec_to_proforma(peptide: str, modifications: str, charge: int | None = None) -> Peptidoform:
    """
    Convert Peptide Record notation to :py:class:`~psm_utils.peptidoform.Peptidoform`.

    Parameters
    ----------
    peptide
        Stripped peptide sequence.
    modifications
        Modifications in Peptide Record notation (e.g., ``4|Oxidation``)
    charge
        Precursor charge state

    Returns
    -------
    peptidoform
        Peptidoform

    Raises
    ------
    InvalidPeprecModificationError
        If a PEPREC modification cannot be parsed.

    """
    # List of peptide sequence with added terminal positions
    peptide_list = [""] + list(peptide) + [""]

    # Add modification labels
    for position, label in zip(modifications.split("|")[::2], modifications.split("|")[1::2]):
        try:
            peptide_list[int(position)] += f"[{label}]"
        except ValueError as e:
            raise InvalidPeprecModificationError(
                f"Could not parse PEPREC modification `{modifications}`."
            ) from e
        except IndexError as e:
            raise InvalidPeprecModificationError(
                f"PEPREC modification has invalid position {position} in peptide "
                f"`{''.join(peptide_list)}`."
            ) from e

    # Add dashes between residues and termini, and join sequence
    peptide_list[0] = peptide_list[0] + "-" if peptide_list[0] else ""
    peptide_list[-1] = "-" + peptide_list[-1] if peptide_list[-1] else ""
    proforma_seq = "".join(peptide_list)

    # Add charge state
    if charge:
        proforma_seq += f"/{charge}"

    return Peptidoform(proforma_seq)


def proforma_to_peprec(peptidoform: Peptidoform) -> tuple[str, str, int | None]:
    """
    Convert :py:class:`~psm_utils.peptidoform.Peptidoform` to Peptide Record notation.

    Parameters
    ----------
    peptidoform
        Input peptidoform object.

    Returns
    -------
    peptide
        Stripped peptide sequence
    modifications
        Modifications in Peptide Record notation
    charge
        Precursor charge state, if available, else :py:const:`None`

    Notes
    -----
    Labile, unlocalized, and fixed modifications are not encoded in the Peptide Record
    notation. To encode fixed modifications, use
    :py:meth:`~psm_utils.psm_list.PSMList.apply_fixed_modifications` before writing to
    Peptide Record.

    """

    def _mod_to_ms2pip(mod_list: list, location: int) -> str:
        """Convert Proforma modification site to MS²PIP modification."""
        if len(mod_list) > 1:
            raise InvalidPeprecModificationError(
                "Multiple modifications per site not supported in Peptide Record format."
            )
        return "|".join([str(location), str(mod_list[0].value)])

    ms2pip_mods = []
    if peptidoform.properties["n_term"]:
        ms2pip_mods.append(_mod_to_ms2pip(peptidoform.properties["n_term"], 0))
    for i, (aa, mod) in enumerate(peptidoform.parsed_sequence):
        if mod:
            ms2pip_mods.append(_mod_to_ms2pip(mod, i + 1))
    if peptidoform.properties["c_term"]:
        ms2pip_mods.append(_mod_to_ms2pip(peptidoform.properties["c_term"], -1))

    peptide = peptidoform.sequence
    modifications = "|".join(ms2pip_mods) if ms2pip_mods else "-"
    if peptidoform.properties["charge_state"]:
        charge = peptidoform.properties["charge_state"].charge
    else:
        charge = None

    return peptide, modifications, charge


def from_dataframe(peprec_df: pd.DataFrame) -> PSMList:
    """
    Convert Peptide Record Pandas DataFrame into PSMList.

    Parameters
    ----------
    peprec_df
        Peptide Record DataFrame

    Returns
    -------
    psm_list
        PSMList object

    """
    psm_list = []
    for _, row in peprec_df.iterrows():
        entry = _PeprecEntry(**row.to_dict())
        psm_list.append(PeptideRecordReader._entry_to_psm(entry, filename=""))
    return PSMList(psm_list=psm_list)


def to_dataframe(psm_list: PSMList) -> pd.DataFrame:
    """
    Convert PSMList object into Peptide Record Pandas DataFrame.

    Parameters
    ----------
    psm_list
        Input PSMList object.

    Returns
    -------
    pd.DataFrame
        Peptide Record DataFrame.

    Examples
    --------
    >>> psm_list = PeptideRecordReader("peprec.csv").read_file()
    >>> psm_utils.io.peptide_record.to_dataframe(psm_list)
        spec_id    peptide               modifications  charge  label  ...
    0  peptide1      ACDEK                           -       2      1  ...
    1  peptide2    ACDEFGR           2|Carbamidomethyl       3      1  ...
    2  peptide3  ACDEFGHIK  0|Acetyl|2|Carbamidomethyl       2      1  ...

    """
    return pd.DataFrame([PeptideRecordWriter._psm_to_entry(psm) for psm in psm_list])


class InvalidPeprecError(PSMUtilsIOException):
    """Invalid Peptide Record file."""


class InvalidPeprecModificationError(InvalidPeprecError):
    """Invalid Peptide Record modification."""
