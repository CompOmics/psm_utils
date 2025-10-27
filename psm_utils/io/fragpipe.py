"""
Reader for PSM files from the Fragpipe platform.

Reads the Philosopher ``psm.tsv`` file as defined on the
`Fragpipe documentation page <https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html>`_.

Notes
-----
- Decoy PSMs and q-values are not returned by FragPipe.

"""

from __future__ import annotations

import csv
from abc import ABC
from collections.abc import Iterator
from pathlib import Path
from typing import Any, cast

import pandas as pd
from pyteomics.proforma import MassModification, to_proforma  # type: ignore[import]

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.io.exceptions import InvalidModificationError
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()

RESCORING_FEATURES: list[str] = [
    "Peptide Length",
    "Retention",
    "Observed Mass",
    "Observed M/Z",
    "Calculated Peptide Mass",
    "Calculated M/Z",
    "Delta Mass",
    "Number of Missed Cleavages",
]


class FragPipeReader(ReaderBase, ABC):
    """Reader for FragPipe TSV format."""

    use_calibrated_mz: bool
    _mz_key: str

    def __init__(
        self,
        filename: str | Path,
        use_calibrated_mz: bool = True,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """
        Reader for MSFragger ``psm.tsv`` file.

        Parameters
        ----------
        filename
            Path to PSM file.
        use_calibrated_mz
            Whether to use ``Calibrated Observed M/Z`` (true) or non-calibrated ``Observed m/z``
            (false), by default True.
        *args
            Additional positional arguments passed to the base class.
        **kwargs
            Additional keyword arguments passed to the base class.

        """
        super().__init__(filename, *args, **kwargs)
        self.use_calibrated_mz = use_calibrated_mz

        self._mz_key = "Calibrated Observed M/Z" if use_calibrated_mz else "Observed M/Z"

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")
            for row in reader:
                yield self._get_peptide_spectrum_match(row)

    def _get_peptide_spectrum_match(self, psm_dict: dict[str, Any]) -> PSM:
        """Parse a single PSM from a FragPipe PSM file."""
        rescoring_features = {ft: psm_dict[ft] for ft in RESCORING_FEATURES if ft in psm_dict}

        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["Peptide"], psm_dict["Assigned Modifications"], psm_dict["Charge"]
            ),
            spectrum_id=self._parse_spectrum_id(psm_dict["Spectrum"]),
            run=self._parse_run(psm_dict["Spectrum File"]),
            is_decoy=False,
            # Assuming this is 1 - PEP, as described in the PeptideProphet paper
            # (https://doi.org/10.1186/1471-2105-13-S16-S1)
            pep=1 - float(psm_dict["Probability"]),
            score=psm_dict["Hyperscore"],
            precursor_mz=psm_dict[self._mz_key],
            retention_time=float(psm_dict["Retention"]),
            ion_mobility=float(psm_dict["Ion Mobility"]) if "Ion Mobility" in psm_dict else None,
            protein_list=self._parse_protein_list(
                psm_dict["Protein"], psm_dict["Mapped Proteins"]
            ),
            source="FragPipe",
            provenance_data=({"fragpipe_filename": str(self.filename)}),
            rescoring_features=rescoring_features,
            metadata={},
        )

    @staticmethod
    def _parse_peptidoform(peptide: str, modifications: str, charge: str | None) -> str:
        """Parse the peptidoform from the modified peptide, peptide, and charge columns."""
        sequence: list[tuple[str, list[MassModification]]] = [(aa, []) for aa in peptide]
        n_term: list[MassModification] = []
        c_term: list[MassModification] = []

        if not modifications:
            return to_proforma(sequence, n_term=n_term, c_term=c_term, charge_state=charge)

        for mod_entry in modifications.split(", "):
            if not mod_entry:
                continue

            parsed_mod_entry: list[str] = mod_entry[:-1].split("(")
            if not len(parsed_mod_entry) == 2:
                raise InvalidModificationError(
                    f"Could not parse modification entry '{mod_entry}'."
                )
            site: str = parsed_mod_entry[0]
            mass: float = float(parsed_mod_entry[1])

            if site == "N-term":
                n_term.append(MassModification(mass))
            elif site == "C-term":
                c_term.append(MassModification(mass))
            else:
                residue: str = site[-1]
                idx: int = int(site[:-1]) - 1
                if idx < 0 or idx >= len(sequence):
                    raise InvalidModificationError(
                        f"Modification position {idx + 1} is out of bounds for peptide of "
                        f"length {len(sequence)}."
                    )
                if sequence[idx][0] != residue:
                    raise InvalidModificationError(
                        f"Modification site residue '{residue}' does not match "
                        f"peptide sequence residue '{sequence[idx][0]}' at position {idx + 1}."
                    )
                sequence[idx][1].append(MassModification(mass))

        return to_proforma(sequence, n_term=n_term, c_term=c_term, charge_state=charge)

    @staticmethod
    def _parse_spectrum_id(spectrum: str) -> str:
        """Extract scan number from spectrum ID: ``(file name).(scan #).(scan #).(charge).``."""
        try:
            return spectrum.split(".")[-2]
        except IndexError:
            return spectrum

    @staticmethod
    def _parse_protein_list(razor_protein: str, mapped_proteins: str | None) -> list[str]:
        """Combine razor protein and mapped proteins into a single list."""
        if mapped_proteins:
            mapped_proteins_list = mapped_proteins.split(", ")
            return [razor_protein] + mapped_proteins_list
        else:
            return [razor_protein]

    @staticmethod
    def _parse_run(spectrum_file: str) -> str:
        """Extract run name from spectrum file."""
        # Depending on the FragPipe workflow used, the run name can be different. In most cases
        # something like 'interact-<run_name>.pep.xml' is used
        if (spectrum_file.endswith(".pep.xml")) and (spectrum_file.startswith("interact-")):
            spectrum_file = spectrum_file.replace("interact-", "")
            return Path(Path(spectrum_file).stem).stem
        else:
            return Path(spectrum_file).stem

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame) -> PSMList:
        """Create a PSMList from a Pandas DataFrame."""
        # Create a temporary reader instance to access the parsing method
        temp_reader = cls(filename="")

        return PSMList(
            psm_list=[
                temp_reader._get_peptide_spectrum_match(cast(dict[str, Any], entry))
                for entry in dataframe.to_dict(orient="records")
            ]
        )
