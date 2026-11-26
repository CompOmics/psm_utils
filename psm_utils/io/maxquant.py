"""Interface to MaxQuant msms.txt PSM files."""

from __future__ import annotations

import csv
import logging
import re
from collections.abc import Callable, Iterator, Mapping, Sequence
from itertools import compress
from pathlib import Path

import numpy as np
from pyteomics import proforma  # type: ignore[import]

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

_NAMED_MODIFICATIONS = {
    "Acetyl (K)": "Acetyl",
    "Acetyl (N-term)": "Acetyl",
    "Acetyl (Peptide N-term)": "Acetyl",
    "Acetyl (Protein N-term)": "Acetyl",
    "Amidated (C-term)": "Amidated",
    "Amidated (Peptide C-term)": "Amidated",
    "Amidated (Protein C-term)": "Amidated",
    "Carbamidomethyl (C)": "Carbamidomethyl",
    "Carbamyl (N-term)": "Carbamyl",
    "Deamidated (NQ)": "Deamidated",
    "Deamidation (N)": "Deamidated",
    "Deamidation (NQ)": "Deamidated",
    "Deamidation (Q)": "Deamidated",
    "Deamidation": "Deamidated",
    "Gln->pyro-Glu": "Gln->pyro-Glu",
    "Glu->pyro-Glu": "Glu->pyro-Glu",
    "GlyGly (K)": "GG",
    "Oxidation (M)": "Oxidation",
    "Oxidation (P)": "Oxidation",
    "Phospho (S)": "Phospho",
    "Phospho (ST)": "Phospho",
    "Phospho (STY)": "Phospho",
    "Phospho (T)": "Phospho",
    "Phospho (Y)": "Phospho",
}

# An underscore denotes an N-terminal tag. Some short codes depend on their site.
_SHORT_MODIFICATIONS = {
    ("ac", "_"): "Acetyl",
    ("ac", "K"): "Acetyl",
    ("ca", "_"): "Carbamyl",
    ("ca", "C"): "Carbamidomethyl",
    ("cm", "C"): "Carbamidomethyl",
    ("de", "N"): "Deamidated",
    ("de", "Q"): "Deamidated",
    ("gl", "_Q"): "Gln->pyro-Glu",
    ("gl", "K"): "GG",
    ("ox", "M"): "Oxidation",
    ("ox", "P"): "Oxidation",
    ("ph", "S"): "Phospho",
    ("ph", "T"): "Phospho",
    ("ph", "Y"): "Phospho",
}


class MSMSReader(ReaderBase):
    """Reader for MaxQuant msms.txt PSM files."""

    def __init__(
        self,
        filename: str | Path,
        *args,
        modification_mapping: Mapping[str, str] | None = None,
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
        modification_mapping
            Optional exact MaxQuant label to ProForma tag mappings. Supplied entries
            override built-in mappings. Pass an empty mapping to keep raw labels.
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
        self._map_modifications = modification_mapping is None or bool(modification_mapping)
        self._modification_overrides = dict(modification_mapping or {})
        self._checked_labels: set[str] = set()
        self._validate_msms()

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        if self._map_modifications:
            logger.info("Applying MaxQuant modification mapping to %s", self.filename)
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
        resolve_label = None
        if self._map_modifications:
            summary = psm_dict.get("Modifications") or ""
            reported = {
                re.sub(r"^\d+\s+", "", item.strip())
                for item in re.split(r"[,;]", summary)
                if item.strip()
            }

            def resolve_label(label: str, site: str) -> str:
                return self._resolve_modification(label, site, reported)

        psm = PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["Modified sequence"], int(psm_dict["Charge"]), resolve_label
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

    def _resolve_modification(self, label: str, site: str, reported: set[str]) -> str:
        """Resolve a literal label using the site and optional row summary."""
        if label in self._modification_overrides:
            return self._modification_overrides[label]
        if label in _NAMED_MODIFICATIONS:
            return _NAMED_MODIFICATIONS[label]

        code = "ox" if label == "Ox" else label
        target = _SHORT_MODIFICATIONS.get((code, site))
        if target is None and site.startswith("_"):
            target = _SHORT_MODIFICATIONS.get((code, "_"))
        if target is not None:
            # Some exports use gl on lysine for a different modification.
            if not reported and (code, site) != ("gl", "K"):
                return target
            for description in reported:
                if _NAMED_MODIFICATIONS.get(description, description) != target:
                    continue
                # Check any explicit site on a recognized description.
                _, _, specificity = description.removesuffix(")").rpartition(" (")
                if specificity.isalpha() and specificity.isupper() and site not in specificity:
                    continue
                if specificity.endswith("N-term") and not site.startswith("_"):
                    continue
                return target

        if label not in self._checked_labels:
            self._checked_labels.add(label)
            try:
                tag = proforma.process_tag_tokens(label)
                if isinstance(tag, proforma.GenericModification):
                    try:
                        proforma.UnimodModification.resolver(
                            name=tag.value, strict=True, exhaustive=False
                        )
                    except KeyError:
                        proforma.PSIModModification.resolver(name=tag.value)
                elif isinstance(tag, proforma.UnimodModification):
                    name, accession = tag.resolver.parse_identifier(tag.value)
                    tag.resolver(name=name, id=accession, strict=True, exhaustive=False)
            except (KeyError, ValueError, OSError, proforma.ProFormaError):
                logger.warning(
                    "Unmapped MaxQuant modification %r in %s. Keeping the original label; "
                    "use modification_mapping to provide a ProForma tag.",
                    label,
                    self.filename,
                )
        return label

    @staticmethod
    def _parse_peptidoform(
        modified_seq: str,
        charge: int,
        resolve_label: Callable[[str, str], str] | None = None,
    ) -> Peptidoform:
        """Parse modified sequence to Peptidoform."""
        modified_seq = modified_seq.strip("_")
        if "(" not in modified_seq:
            return Peptidoform(f"{modified_seq}/{charge}")

        # Split the unmodified spans and literal labels without interpreting nested
        # parentheses or using label text as part of a regular expression.
        segments: list[str] = []
        labels: list[str] = []
        depth = start = previous = 0
        for index, char in enumerate(modified_seq):
            if char == "(":
                if depth == 0:
                    start = index
                    segments.append(modified_seq[previous:start])
                depth += 1
            elif char == ")":
                depth -= 1
                if depth == 0:
                    labels.append(modified_seq[start + 1 : index])
                    previous = index + 1
                elif depth < 0:
                    break
        if depth:
            raise MSMSParsingError(f"Unbalanced modification in {modified_seq!r}.")
        segments.append(modified_seq[previous:])

        site = "_" + next((segment[0] for segment in segments if segment), "")
        last_segment = max((i for i, segment in enumerate(segments) if segment), default=0)
        sequence, n_term, c_term = [], [], []
        for index, (segment, label) in enumerate(zip(segments, labels)):
            sequence.append(segment)
            if segment:
                site = segment[-1]
            tag = f"[{resolve_label(label, site) if resolve_label is not None else label}]"
            if site.startswith("_"):
                n_term.append(tag)
            elif index >= last_segment and (
                label == "Amidated"
                or label.endswith((" (C-term)", " (Peptide C-term)", " (Protein C-term)"))
                or (resolve_label is None and index == len(labels) - 1)
            ):
                c_term.append(tag)
            else:
                sequence.append(tag)

        if n_term:
            sequence[:0] = [*n_term, "-"]
        sequence.append(segments[-1])
        if c_term:
            sequence.extend(["-", *c_term])
        return Peptidoform(f"{''.join(sequence)}/{charge}")


class MSMSParsingError(PSMUtilsException):
    """Error while parsing MaxQuant msms.txt PSM file."""
