"""
Interface with OpenMS idXML PSM files.

Notes
-----
* idXML supports multiple peptide hits (identifications) per spectrum. Each peptide hit
  is parsed as an individual :py:class:`~psm_utils.psm.PSM` object.

"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any, cast
from warnings import filterwarnings

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

filterwarnings(
    "ignore",
    message="OPENMS_DATA_PATH environment variable already exists",
    category=UserWarning,
    module="pyopenms",
)

try:
    import pyopenms as oms  # type: ignore[import]

    _has_openms = True
except ImportError:
    _has_openms = False
    oms = None  # type: ignore[assignment]

logger = logging.getLogger(__name__)

DEFAULT_SCORE_TYPE = "search_engine_score"
TARGET_DECOY_KEY = "target_decoy"
QVALUE_KEY = "q-value"
PEP_KEY = "PEP"
SPECTRUM_REFERENCE_KEY = "spectrum_reference"
ID_MERGE_INDEX_KEY = "id_merge_index"
SPECTRA_DATA_KEY = "spectra_data"
ION_MOBILITY_KEY = "IM"

# Patterns to match open and closed round/square brackets
MOD_PATTERN = re.compile(r"\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)")
MOD_PATTERN_NTERM = re.compile(r"^\.\[((?:[^][]+|\[(?:[^][]+|\[[^][]*\])*\])*)\]")
MOD_PATTERN_CTERM = re.compile(r"\.\[((?:[^][]+|\[(?:[^][]+|\[[^][]*\])*\])*)\]$")

# Extracted from the OpenMS PSMFeatureExtractor, which adds and manipulates features that will be given to percolator
# https://github.com/OpenMS/OpenMS/blob/342f6524e76a2bab3dcb428ba2f4aa2d6bfe8483/src/topp/PSMFeatureExtractor.cpp
RESCORING_FEATURE_LIST = [
    "isotope_error",
    "MS:1002049",  # MSGFPlus unchanged RawScore
    "MS:1002050",  # MSGFPlus unchanged DeNovoScore
    "MSGF:ScoreRatio",
    "MSGF:Energy",
    "MSGF:lnEValue",
    "MSGF:lnExplainedIonCurrentRatio",
    "MSGF:lnNTermIonCurrentRatio",
    "MSGF:lnCTermIonCurrentRatio",
    "MSGF:lnMS2IonCurrent",
    "MSGF:MeanErrorTop7",
    "MSGF:sqMeanErrorTop7",
    "MSGF:StdevErrorTop7",
    "XTANDEM:hyperscore",
    "XTANDEM:deltascore",
    "MS:1001330",  # expect_score
    "hyperscore",  # MSFragger
    "nextscore",  # MSFragger
    "COMET:deltaCn",  # recalculated deltaCn = (current_XCorr - 2nd_best_XCorr) / max(current_XCorr, 1)
    "COMET:deltaLCn",  # deltaLCn = (current_XCorr - worst_XCorr) / max(current_XCorr, 1)
    "COMET:lnExpect",  # log(E-value)
    "MS:1002252",  # unchanged XCorr
    "MS:1002255",  # unchanged Sp = number of candidate peptides
    "COMET:lnNumSP",  # log(number of candidate peptides)
    "COMET:lnRankSP",  # log(rank based on Sp score)
    "COMET:IonFrac",  # matched_ions / total_ions
    "MS:1001171",  # unchanged mScore
    "MASCOT:delta_score",  # delta score based on mScore
    "CONCAT:lnEvalue",
    "CONCAT:deltaLnEvalue",
    "SAGE:ln(-poisson)",
    "SAGE:ln(delta_best)",
    "SAGE:ln(delta_next)",
    "SAGE:ln(matched_intensity_pct)",
    "SAGE:longest_b",
    "SAGE:longest_y",
    "SAGE:longest_y_pct",
    "SAGE:matched_peaks",
    "SAGE:scored_candidates",
]


class IdXMLReader(ReaderBase):
    """Reader for idXML files with comprehensive type safety and error handling."""

    protein_ids: Any  # list[oms.ProteinIdentification]
    peptide_ids: Any  # list[oms.PeptideIdentification]
    user_params_metadata: list[str]
    rescoring_features: list[str]

    def __init__(self, filename: Path | str, *args: Any, **kwargs: Any) -> None:
        """
        Reader for idXML files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to idXML file.
        *args
            Additional positional arguments passed to the base class.
        **kwargs
            Additional keyword arguments passed to the base class.

        Examples
        --------
        >>> from psm_utils.io import IdXMLReader
        >>> reader = IdXMLReader("example.idXML")
        >>> psm_list = [psm for psm in reader]

        """
        super().__init__(filename, *args, **kwargs)
        if not _has_openms:
            raise ImportError("pyOpenMS is required to read idXML files")

        self.protein_ids, self.peptide_ids = self._parse_idxml()
        self.user_params_metadata = self._get_userparams_metadata(self.peptide_ids[0].getHits()[0])
        self.rescoring_features = self._get_rescoring_features(self.peptide_ids[0].getHits()[0])

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        for peptide_id in self.peptide_ids:
            for peptide_hit in peptide_id.getHits():
                yield self._parse_psm(self.protein_ids, peptide_id, peptide_hit)

    def _parse_idxml(self) -> tuple[Any, Any]:
        """
        Parse idXML using pyopenms and perform sanity checks to make sure the file is not empty.

        Returns
        -------
        tuple of (Any, Any)
            Tuple containing (ProteinIdentification, PeptideIdentification) lists

        Raises
        ------
        IdXMLReaderEmptyListException
            If the idXML file contains no data to parse

        """
        protein_ids: Any = []  # list[oms.ProteinIdentification]
        peptide_ids: Any = []  # list[oms.PeptideIdentification]
        oms.IdXMLFile().load(str(self.filename), protein_ids, peptide_ids)  # type: ignore

        if len(protein_ids) == 0:
            raise IdXMLReaderEmptyListException(
                f"File {self.filename} contains no proteins. Nothing to parse."
            )
        elif len(peptide_ids) == 0:
            raise IdXMLReaderEmptyListException(
                f"File {self.filename} contains no PeptideIdentifications. Nothing to parse."
            )
        elif len(peptide_ids[0].getHits()) == 0:
            raise IdXMLReaderEmptyListException(
                f"File {self.filename} contains no PeptideHits. Nothing to parse."
            )
        else:
            return protein_ids, peptide_ids

    @staticmethod
    def _parse_peptidoform(sequence: str, charge: int) -> str:
        """
        Parse idXML peptide to :py:class:`~psm_utils.peptidoform.Peptidoform`.

        Parameters
        ----------
        sequence
            Peptide sequence in idXML format
        charge
            Precursor charge state

        Returns
        -------
        str
            Peptide sequence in Peptidoform format with charge

        Notes
        -----
        Implemented according to the documentation on
        `github.com/OpenMS/OpenMS <https://github.com/OpenMS/OpenMS/blob/8cb90/src/openms/include/OpenMS/CHEMISTRY/AASequence.h>`_
        . The differentiation between square- and round bracket notation is removed after parsing.

        """
        sequence = MOD_PATTERN.sub(r"[\1]", sequence)
        if sequence[:2] == ".[":
            sequence = MOD_PATTERN_NTERM.sub(r"[\1]-", sequence)
        if sequence[-1] == "]":
            sequence = MOD_PATTERN_CTERM.sub(r"-[\1]", sequence)
        sequence = sequence.strip(".")
        sequence += f"/{charge}"

        return sequence

    def _parse_psm(
        self,
        protein_ids: Any,
        peptide_id: Any,
        peptide_hit: Any,
    ) -> PSM:
        """
        Parse idXML :py:class:`~pyopenms.PeptideHit` to :py:class:`~psm_utils.psm.PSM`.

        Uses additional information from :py:class:`~pyopenms.ProteinIdentification` and
        :py:class:`~pyopenms.PeptideIdentification` to annotate parameters of the
        :py:class:`~psm_utils.psm.PSM` object.

        Parameters
        ----------
        protein_ids
            List of ProteinIdentification objects
        peptide_id
            PeptideIdentification object
        peptide_hit
            PeptideHit object

        Returns
        -------
        PSM
            Parsed PSM object with all available information

        """
        peptidoform = self._parse_peptidoform(
            peptide_hit.getSequence().toString(), peptide_hit.getCharge()
        )
        # This is needed to calculate a qvalue before rescoring the PSMList
        peptide_id_metadata = {
            "idxml:score_type": str(peptide_id.getScoreType()),
            "idxml:higher_score_better": str(peptide_id.isHigherScoreBetter()),
            "idxml:significance_threshold": str(peptide_id.getSignificanceThreshold()),
        }
        peptide_hit_metadata = {
            key: peptide_hit.getMetaValue(key) for key in self.user_params_metadata
        }
        return PSM(
            peptidoform=peptidoform,
            spectrum_id=peptide_id.getMetaValue(SPECTRUM_REFERENCE_KEY),
            run=self._get_run(protein_ids, peptide_id),
            is_decoy=self._is_decoy(peptide_hit),
            score=peptide_hit.getScore(),
            precursor_mz=peptide_id.getMZ(),
            retention_time=peptide_id.getRT(),
            ion_mobility=self._get_ion_mobility(peptide_hit),
            protein_list=[
                accession.decode() for accession in peptide_hit.extractProteinAccessionsSet()
            ],
            rank=peptide_hit.getRank() + 1,  # 0-based to 1-based
            source="idXML",
            # Storing proforma notation of peptidoform and UNIMOD peptide sequence for mapping back
            # to original sequence in writer
            provenance_data={str(peptidoform): peptide_hit.getSequence().toString()},
            # Store metadata of PeptideIdentification and PeptideHit objects
            metadata={**peptide_id_metadata, **peptide_hit_metadata},
            rescoring_features={
                key: float(peptide_hit.getMetaValue(key))  # type: ignore
                for key in self.rescoring_features
            },
        )

    @staticmethod
    def _get_run(protein_ids: Any, peptide_id: Any) -> str | None:
        """
        Get run name from idXML using pyopenms.

        If the idXML file contains a merge index, use it to annotate the run name without file
        extension.
        """
        # Check if spectra_data is available
        if not protein_ids[0].metaValueExists(SPECTRA_DATA_KEY):
            return None

        spectra_data = cast(list[bytes], protein_ids[0].getMetaValue(SPECTRA_DATA_KEY))

        # Determine index to use
        if peptide_id.metaValueExists(ID_MERGE_INDEX_KEY):
            index = cast(int, peptide_id.getMetaValue(ID_MERGE_INDEX_KEY))
        else:
            index = 0

        # Extract run path
        try:
            run_path = Path(spectra_data[index].decode()).stem
        except (IndexError, UnicodeDecodeError):
            return None

        # Handle the special case where run path is the string "None"
        return None if run_path == "None" else run_path

    @staticmethod
    def _get_ion_mobility(peptide_hit: Any) -> float | None:
        """
        Get ion mobility from PeptideHit.

        Parameters
        ----------
        peptide_hit
            PeptideHit object

        Returns
        -------
        float or None
            Ion mobility value or None if not available or invalid

        """
        if not peptide_hit.metaValueExists(ION_MOBILITY_KEY):
            return None

        im_value = peptide_hit.getMetaValue(ION_MOBILITY_KEY)
        try:
            return float(im_value)  # type: ignore[arg-type]
        except (ValueError, TypeError):
            return None

    def _get_userparams_metadata(self, peptide_hit: Any) -> list[str]:
        """Get list of string type UserParams attached to each PeptideHit."""
        # Fill the key list with all the keys from the PeptideHit
        # Empty list is required for the Cython wrapper to work correctly
        keys: list[bytes] = []
        peptide_hit.getKeys(keys)

        return [
            key.decode()
            for key in keys
            if not self._is_float(peptide_hit.getMetaValue(key.decode()))
        ]

    def _get_rescoring_features(self, peptide_hit: Any) -> list[str]:
        """Get list of rescoring features in UserParams attached to each PeptideHit."""
        keys: list[bytes] = []
        peptide_hit.getKeys(keys)

        return [
            key.decode()
            for key in keys
            if (
                self._is_float(peptide_hit.getMetaValue(key.decode()))
                and key.decode() in RESCORING_FEATURE_LIST
            )
        ]

    @staticmethod
    def _is_float(element: Any) -> bool:
        """Check if element can be coerced to a float."""
        if element is None:
            return False
        try:
            float(element)
            return True
        except ValueError:
            return False

    @staticmethod
    def _is_decoy(peptide_hit: Any) -> bool | None:
        """Check if PSM is target or decoy."""
        if peptide_hit.metaValueExists(TARGET_DECOY_KEY):
            return peptide_hit.getMetaValue(TARGET_DECOY_KEY) == "decoy"
        else:
            return None


class IdXMLWriter(WriterBase):
    """Writer for idXML files with comprehensive error handling."""

    protein_ids: Any | None
    peptide_ids: Any | None

    def __init__(
        self,
        filename: str | Path,
        *args: Any,
        protein_ids: Any | None = None,
        peptide_ids: Any | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Writer for idXML files.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments passed to the base class.
        protein_ids
            Optional list of :py:class:`~pyopenms.ProteinIdentification` objects to be written to the idXML file.
        peptide_ids
            Optional list of :py:class:`~pyopenms.PeptideIdentification` objects to be written to the idXML file.
        **kwargs
            Additional keyword arguments passed to the base class.

        Notes
        -----
        - Unlike other psm_utils.io writer classes, :py:class:`IdXMLWriter` does not support writing
          a single PSM to a file with the :py:meth:`write_psm` method. Only writing a full PSMList
          to a file at once with the :py:meth:`write_file` method is currently supported.
        - If `protein_ids` and `peptide_ids` are provided, each :py:class:`~pyopenms.PeptideIdentification`
          object in the list `peptide_ids` will be updated with new `rescoring_features` from the PSMList.
          Otherwise, new pyopenms objects will be created, filled with information of PSMList and written to the idXML file.

        Examples
        --------
        - Example with `pyopenms` objects:

        >>> from psm_utils.io.idxml import IdXMLReader, IdXMLWriter
        >>> reader = IdXMLReader("psm_utils/tests/test_data/test_in.idXML")
        >>> psm_list = reader.read_file()
        >>> for psm in psm_list:
        ...     psm.rescoring_features = {**psm.rescoring_features, **{"feature": 1}}
        >>> writer = IdXMLWriter("psm_utils/tests/test_data//test_out.idXML", reader.protein_ids, reader.peptide_ids)
        >>> writer.write_file(psm_list)

        - Example without `pyopenms` objects:

        >>> from psm_utils.psm_list import PSMList
        >>> psm_list = PSMList(psm_list=[PSM(peptidoform="ACDK", spectrum_id=1, score=140.2, retention_time=600.2)])
        >>> writer = IdXMLWriter("psm_utils/tests/test_data//test_out.idXML")
        >>> writer.write_file(psm_list)

        """
        super().__init__(filename, *args, **kwargs)
        if not _has_openms:
            raise ImportError("pyOpenMS is required to write idXML files")

        self.protein_ids = protein_ids
        self.peptide_ids = peptide_ids

    def __enter__(self) -> IdXMLWriter:
        """Open file for writing and return self."""
        return self

    def __exit__(self, *args, **kwargs) -> None:
        """Close file and writer."""
        pass

    def write_psm(self, psm: PSM):
        """
        Write a single PSM to the PSM file.

        This method is currently not supported (see Notes).

        Parameters
        ----------
        psm
            PSM object to write

        Raises
        ------
        NotImplementedError
            IdXMLWriter currently does not support write_psm.

        """
        raise NotImplementedError("IdXMLWriter currently does not support write_psm.")

    def write_file(self, psm_list: PSMList) -> None:
        """
        Write the PSMList to the PSM file.

        If `self.protein_ids` and `self.peptide_ids` are not None, the PSM list scores, ranks,
        and rescoring features will first be merged with the existing IDs from those objects.

        Parameters
        ----------
        psm_list
            List of PSM objects to write to file

        """
        psm_dict = psm_list.get_psm_dict()

        if self.protein_ids is not None and self.peptide_ids is not None:
            self._update_existing_ids(psm_dict)
        # Check if one of self.protein_ids or self.peptide_ids is None
        elif self.protein_ids is not None or self.peptide_ids is not None:
            logger.warning(
                "One of the protein_ids or peptide_ids is None. Falling back to creating new "
                "idXML files solely based on the PSMList."
            )
            self._create_new_ids(psm_dict)
        else:
            self._create_new_ids(psm_dict)

    def _update_existing_ids(
        self, psm_dict: dict[str | None, dict[str, dict[str, list[PSM]]]]
    ) -> None:
        """
        Update an existing idXML file with info from the PSM list or write a new one.

        Update existing :py:class:`~pyopenms.ProteinIdentification` and
        :py:class:`~pyopenms.PeptideIdentification` objects with new features from the PSMList
        or create new ones.
        """
        if not self.protein_ids or not self.peptide_ids:
            raise IdXMLException(
                "Both protein_ids and peptide_ids must be provided to update existing idXML."
            )
        # Access run name(s) from ProteinIdentification
        spectrum_files = [
            Path(run.decode()).stem
            for run in cast(list[bytes], self.protein_ids[0].getMetaValue(SPECTRA_DATA_KEY))
        ]
        for peptide_id in self.peptide_ids:
            if len(spectrum_files) > 1:
                id_merge_index = cast(int, peptide_id.getMetaValue(ID_MERGE_INDEX_KEY))
                run = spectrum_files[id_merge_index]
            else:
                run = spectrum_files[0]

            # Get PSM objects associated from runs since we are writing a merged idXML
            # NOTE: Collections with multiple protein_ids and peptide_ids is not supported
            try:
                spectrum_ref = cast(str, peptide_id.getMetaValue(SPECTRUM_REFERENCE_KEY))
                psms = psm_dict[None][run][spectrum_ref]
            except KeyError as e:
                raise IdXMLException(
                    "Multiple collections are not supported when parsing single pyopenms protein "
                    "and peptide objects."
                ) from e

            # Dict of UNIMOD peptide sequence and PSM object
            hit_dict = {
                (psm.provenance_data or {})[str(psm.peptidoform)]: psm
                for psm in psms
                if psm.provenance_data and str(psm.peptidoform) in psm.provenance_data
            }
            # Update PeptideHits according to the PSM objects
            updated_peptide_hits = []
            for peptide_hit in peptide_id.getHits():
                sequence = peptide_hit.getSequence().toString()
                psm = hit_dict[sequence]
                self._update_peptide_hit(peptide_hit, psm)
                updated_peptide_hits.append(peptide_hit)

            peptide_id.setHits(updated_peptide_hits)

        oms.IdXMLFile().store(str(self.filename), self.protein_ids, self.peptide_ids)  # type: ignore

    def _update_peptide_hit(self, peptide_hit: Any, psm: PSM) -> None:
        """Inplace update of PeptideHit with novel predicted features information from PSM."""
        # Update core PSM attributes
        if psm.score is not None:
            peptide_hit.setScore(psm.score)
        if psm.rank is not None:
            peptide_hit.setRank(psm.rank - 1)  # 1-based to 0-based
        if psm.qvalue is not None:
            peptide_hit.setMetaValue(QVALUE_KEY, psm.qvalue)
        if psm.pep is not None:
            peptide_hit.setMetaValue(PEP_KEY, psm.pep)

        # Add rescoring features (only those not in the standard list)
        if psm.rescoring_features:
            for feature, value in psm.rescoring_features.items():
                if feature not in RESCORING_FEATURE_LIST:
                    # Convert numpy objects to floats as pyopenms does not support numpy objects
                    peptide_hit.setMetaValue(feature, float(value))

    def _create_new_ids(self, psm_dict: dict[str | None, dict[str, dict[str, list[PSM]]]]) -> None:
        """Create new ProteinIdentification and PeptideIdentification objects with new features."""
        for collection, runs in psm_dict.items():
            self._create_ids_for_collection(collection, runs)

    def _create_ids_for_collection(
        self, collection: str | None, runs: dict[str, dict[str, list[PSM]]]
    ) -> None:
        """Create ProteinIdentification and PeptideIdentification objects for a single collection."""
        self.protein_ids = [oms.ProteinIdentification()]  # type: ignore
        self.peptide_ids = []

        # Set msrun filename with spectra_data meta value
        msrun_reference = [str(run).encode() for run in runs.keys()]
        self.protein_ids[0].setMetaValue(SPECTRA_DATA_KEY, msrun_reference)

        protein_list: list[list[str]] = []

        for run, psm_dict_run in runs.items():
            for spectrum_id, psms in psm_dict_run.items():
                # Collect protein accessions
                protein_list.append(
                    [accession for psm in psms for accession in (psm.protein_list or [])]
                )

                # Create PeptideIdentification
                peptide_id = self._create_peptide_identification(
                    spectrum_id, run, msrun_reference, psms
                )

                # Create PeptideHits
                peptide_hits = [self._create_peptide_hit(psm) for psm in psms]
                peptide_id.setHits(peptide_hits)
                self.peptide_ids.append(peptide_id)

        # Create protein hits
        self._create_protein_hits(protein_list)

        # Write idXML file
        filename = "/".join(filter(None, [collection, str(self.filename)]))
        oms.IdXMLFile().store(filename, self.protein_ids, self.peptide_ids)  # type: ignore

    def _create_peptide_identification(
        self,
        spectrum_id: str,
        run: str,
        msrun_reference: list[bytes],
        psms: list[PSM],
    ) -> Any:
        """Create a PeptideIdentification object for a spectrum."""
        peptide_id = oms.PeptideIdentification()  # type: ignore
        peptide_id.setMetaValue(SPECTRUM_REFERENCE_KEY, spectrum_id)
        peptide_id.setMetaValue(ID_MERGE_INDEX_KEY, msrun_reference.index(str(run).encode()))

        # Set properties from first PSM
        first_psm = psms[0]
        if first_psm.score is not None:
            peptide_id.setScoreType(DEFAULT_SCORE_TYPE)
        if first_psm.precursor_mz is not None:
            peptide_id.setMZ(first_psm.precursor_mz)
        if first_psm.retention_time is not None:
            peptide_id.setRT(first_psm.retention_time)

        return peptide_id

    def _create_peptide_hit(self, psm: PSM) -> Any:
        """Create a PeptideHit object from a PSM."""
        peptide_hit = oms.PeptideHit()  # type: ignore

        # Set sequence
        peptide_hit.setSequence(
            oms.AASequence.fromString(  # type: ignore
                self._convert_proforma_to_unimod(psm.peptidoform)
            )
        )

        # Set charge
        if psm.peptidoform.precursor_charge is not None:
            peptide_hit.setCharge(psm.peptidoform.precursor_charge)

        # Set target/decoy information
        target_decoy_value = (
            "" if psm.is_decoy is None else ("decoy" if psm.is_decoy else "target")
        )
        peptide_hit.setMetaValue(TARGET_DECOY_KEY, target_decoy_value)

        # Set optional values
        if psm.qvalue is not None:
            peptide_hit.setMetaValue(QVALUE_KEY, psm.qvalue)
        if psm.pep is not None:
            peptide_hit.setMetaValue(PEP_KEY, psm.pep)
        if psm.rank is not None:
            peptide_hit.setRank(psm.rank - 1)  # 1-based to 0-based

        # Add metadata and features
        if psm.metadata:
            self._add_meta_values_from_dict(peptide_hit, psm.metadata)
        if psm.provenance_data:
            self._add_meta_values_from_dict(peptide_hit, psm.provenance_data)
        if psm.rescoring_features:
            self._add_meta_values_from_dict(peptide_hit, psm.rescoring_features)

        # Add protein evidence
        if psm.protein_list is not None:
            for protein in psm.protein_list:
                peptide_evidence = oms.PeptideEvidence()  # type: ignore
                peptide_evidence.setProteinAccession(protein)
                peptide_hit.addPeptideEvidence(peptide_evidence)

        return peptide_hit

    def _create_protein_hits(self, protein_list: list[list[str]]) -> None:
        """Create protein hits from collected protein accessions."""
        # Get unique protein accessions
        unique_proteins = list(
            {accession for protein_sublist in protein_list for accession in protein_sublist}
        )

        protein_hits = []
        for accession in unique_proteins:
            protein_hit = oms.ProteinHit()  # type: ignore
            protein_hit.setAccession(accession)
            protein_hits.append(protein_hit)

        if self.protein_ids and len(self.protein_ids) > 0:
            self.protein_ids[0].setHits(protein_hits)

    def _convert_proforma_to_unimod(self, peptidoform: Peptidoform) -> str:
        """Convert a peptidoform sequence in proforma notation to UNIMOD notation."""
        sequence = str(peptidoform).split("/")[0]

        # Replace square brackets around modifications with parentheses
        sequence = re.sub(r"\[([^\]]+)\]", r"(\1)", sequence)

        # Check for N-terminal and C-terminal modifications
        if sequence.startswith("["):
            sequence = re.sub(r"^\[([^\]]+)\]-", r"(\1)", sequence)
        if sequence.endswith("]"):
            sequence = re.sub(r"-\[([^\]]+)\]$", r"-(\1)", sequence)

        # Remove dashes for N-terminal and C-terminal modifications
        sequence = sequence.replace(")-", ")").replace("-(", "(")

        return sequence

    def _add_meta_values_from_dict(self, peptide_hit: Any, d: dict[str, Any] | None) -> None:
        """Add meta values inplace to :py:class:`~pyopenms.PeptideHit` from a dictionary."""
        if d is None:
            return

        for key, value in d.items():
            # Convert numpy objects to floats since pyopenms does not support numpy objects
            if not isinstance(value, str):
                try:
                    value = float(value)
                except (ValueError, TypeError):
                    # Skip values that cannot be converted
                    continue
            peptide_hit.setMetaValue(key, value)


class IdXMLException(PSMUtilsIOException):
    """Exception in psm_utils.io.IdXML."""

    pass


class IdXMLReaderEmptyListException(PSMUtilsIOException):
    """Exception in psm_utils.io.IdXMLReader."""

    pass
