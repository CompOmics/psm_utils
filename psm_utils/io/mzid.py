"""
Reader and writer for HUPO-PSI mzIdentML format PSM files.

See `psidev.info/mzidentml <https://psidev.info/mzidentml>`_ for more info on the
format.

"""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from pathlib import Path
from typing import Any, cast

from psims.mzid import MzIdentMLWriter  # type: ignore[import]
from pyteomics import mzid  # type: ignore[import]
from rich.progress import Progress

from psm_utils import __version__
from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)

# Excerpt from MS:1001143 items (PSM-level search engine specific statistic)
# Not all child terms are used, as not all statistics are direct scores.
# Items are sorted by priority (if more scores are present, the first found one is used)
STANDARD_SEARCHENGINE_SCORES = [
    "peptideshaker psm score",
    "amanda:amandascore",
    "andromeda:score",
    "byonic:score",
    "comet:xcorr",
    "debunker:score",
    "identitye score",
    "ksdp score",
    "ms-gf:specevalue",
    "ms-gf:evalue",
    "ms-gf:rawscore",
    "ms-gf:denovoscore",
    "msfit:mowse score",
    "mspathfinder:rawscore",
    "mspepsearch:score",
    "mascot:score",
    "metamorpheus:score",
    "omssa:evalue",
    "openpepxl:score",
    "peaks:peptidescore",
    "phenyx:pepzscore",
    "prolucid:xcorr",
    "prosight:specral c-score",
    "profound:z value",
    "proteinprospector:score",
    "proteinscape:sequestmetascore",
    "proteomediscoverer:delta score",
    "proteome discoverer delta score",
    "sequest:xcorr",
    "sim-xl score ",
    "sqid:score ",
    "sonar:score",
    "spectrummill:score",
    "topmg:spectral e-value",
    "x!tandem:hyperscore",
    "zcore:probscore:",
    "percolator:score",
    "xi:score",
    "search engine specific score",
]
Q_VALUE_TERMS = [
    "PSM-level q-value",
    "MS-GF:QValue",
    "MSPathFinder:QValue",
    "ProSight:spectral Q-value",
]
PEP_TERMS = [
    "PSM-level local FDR",
    "Andromeda:PEP",
    "Byonic:PEP",
    "MaxQuant-DIA PEP",
    "MS-GF:PEP",
    "OpenMS:ConsensusID PEP",
    "percolator:PEP",
]


class MzidReader(ReaderBase):
    """Reader for HUPO-PSI mzIdentML format PSM files."""

    def __init__(
        self, filename: str | Path, *args: Any, score_key: str | None = None, **kwargs: Any
    ) -> None:
        """
        Reader for HUPO-PSI mzIdentML format PSM files.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments passed to parent class.
        score_key
            Name of the score metric to use as PSM score. If not provided, the score metric is
            inferred from the file if one of the child parameters of ``MS:1001143`` is present.
        **kwargs
            Additional keyword arguments passed to parent class.

        Examples
        --------
        MzidReader supports iteration:

        >>> from psm_utils.io.mzid import MzidReader
        >>> for psm in MzidReader("peptides_1_1_0.mzid"):
        ...     print(psm.peptidoform.proforma)
        ACDEK
        AC[Carbamidomethyl]DEFGR
        [Acetyl]-AC[Carbamidomethyl]DEFGHIK

        Or a full file can be read at once into a :py:class:`psm_utils.psm_list.PSMList`
        object:

        >>> mzid_reader = MzidReader("peptides_1_1_0.mzid")
        >>> psm_list = mzid_reader.read_file()

        Notes
        -----
        - :py:class:`MzidReader` looks for the ``retention time`` or ``scan start time`` cvParams
          in both SpectrumIdentificationResult and SpectrumIdentificationItem levels. Note that
          according to the mzIdentML specification document (v1.1.1) neither cvParams are expected
          to be present at either levels.
        - For the :py:attr:`PSM.spectrum_id` property, the ``spectrum title`` cvParam is preferred
          over the ``spectrumID`` attribute, as these titles always match the titles in the  peak
          list files. ``spectrumID`` is then saved in ``PSM.metadata["mzid_spectrum_id"]``.
          If ``spectrum title`` is absent, ``spectrumID`` is saved to :py:attr:`PSM.spectrum_id`.

        """
        super().__init__(filename, *args, **kwargs)
        self._non_metadata_keys = ["ContactRole", "passThreshold"]
        self._score_key = score_key
        self._rt_key: str | None = None
        self._spectrum_rt_key: str | None = None
        self._qvalue_key: str | None = None
        self._pep_key: str | None = None
        self._im_key: str | None = None

        self._source = self._infer_source()

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over mzIdentML file and return PSMs one-by-one."""
        with mzid.MzIdentML(str(self.filename)) as reader:
            first_entry = next(reader)
            # Parse spectrum metadata
            self._get_toplevel_non_metadata_keys(first_entry.keys())
            # Parse PSM non-metadata keys, rt key and score key
            self._get_non_metadata_keys(first_entry["SpectrumIdentificationItem"][0].keys())

        with mzid.MzIdentML(str(self.filename)) as reader:
            for spectrum in reader:
                # Parse spectrum metadata
                spectrum_id = spectrum["spectrumID"]
                spectrum_title = spectrum.get("spectrum title")
                run = Path(spectrum["location"]).stem if spectrum.get("location") else None
                rt = (
                    float(spectrum[self._spectrum_rt_key])
                    if self._spectrum_rt_key and self._spectrum_rt_key in spectrum
                    else None
                )
                ion_mobility = (
                    float(spectrum[self._im_key])
                    if self._im_key and self._im_key in spectrum
                    else None
                )

                # Parse PSMs from spectrum
                for entry in spectrum["SpectrumIdentificationItem"]:
                    yield self._get_peptide_spectrum_match(
                        spectrum_id, spectrum_title, run, rt, ion_mobility, entry
                    )

    @staticmethod
    def _get_xml_namespace(root_tag: str) -> str:
        """Extract XML namespace from root tag."""
        m = re.match(r"\{.*\}", root_tag)
        return m.group(0) if m else ""

    def _infer_source(self) -> str | None:
        """Infer search engine source from mzIdentML file metadata."""
        mzid_xml = ET.parse(self.filename)
        root = mzid_xml.getroot()
        name_space = self._get_xml_namespace(root.tag)
        software = root.find(f".//{name_space}AnalysisSoftware")
        if software is None:
            return None
        try:
            return software.attrib["name"]
        except KeyError:
            return None

    @staticmethod
    def _parse_peptidoform(
        seq: str, modification_list: list[dict[str, Any]], charge: int | None
    ) -> Peptidoform:
        """Parse mzIdentML sequence and modifications into Peptidoform object."""
        peptide = [""] + list(seq) + [""]

        # Add modification labels
        for mod in modification_list:
            peptide[int(mod["location"])] += f"[{mod['name']}]"

        # Add dashes between residues and termini, and join sequence
        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""
        proforma_seq = "".join(peptide)

        # Add charge state
        if charge:
            proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)

    @staticmethod
    def _parse_peptide_evidence_ref(
        peptide_evidence_list: list[dict[str, Any]],
    ) -> tuple[bool, list[str]]:
        """
        Parse PeptideEvidence list to determine decoy status and protein accessions.

        PSM is considered decoy only if ALL PeptideEvidence entries are decoy.
        """
        is_decoy = all(entry.get("isDecoy", False) for entry in peptide_evidence_list)
        protein_list = [
            accession
            for d in peptide_evidence_list
            if (accession := d.get("accession")) is not None
        ]
        return is_decoy, protein_list

    def _get_peptide_spectrum_match(
        self,
        spectrum_id: str,
        spectrum_title: str | None,
        run: str | None,
        rt: float | None,
        ion_mobility: float | None,
        spectrum_identification_item: dict[str, Any],
    ) -> PSM:
        """Parse single mzIdentML SpectrumIdentificationItem into PSM object."""
        sii = spectrum_identification_item
        modifications = cast(list[dict[str, Any]], sii.get("Modification", []))
        sequence = cast(str, sii["PeptideSequence"])
        charge_value = sii.get("chargeState")
        charge = int(charge_value) if charge_value is not None else None
        peptidoform = self._parse_peptidoform(sequence, modifications, charge)
        is_decoy, protein_list = self._parse_peptide_evidence_ref(
            cast(list[dict[str, Any]], sii["PeptideEvidenceRef"])
        )

        precursor_mz_value = sii.get("experimentalMassToCharge")
        precursor_mz = float(precursor_mz_value) if precursor_mz_value is not None else None

        # Override spectrum-level RT if present at PSM level
        if self._rt_key and self._rt_key in sii:
            rt = float(sii[self._rt_key])

        metadata = {col: str(sii[col]) for col in sii.keys() if col not in self._non_metadata_keys}

        # Prefer spectrum title (identical to titles in MGF), fall back to spectrumID
        if spectrum_title:
            metadata["mzid_spectrum_id"] = spectrum_id
            psm_spectrum_id = spectrum_title
        else:
            psm_spectrum_id = spectrum_id

        score = None
        if self._score_key:
            score_value = sii.get(self._score_key)
            score = float(score_value) if score_value is not None else None

        # Calculate qvalue and pep with cleaner logic
        qvalue = None
        if self._qvalue_key:
            qvalue_raw = sii.get(self._qvalue_key)
            qvalue = float(qvalue_raw) if qvalue_raw is not None else None

        pep = None
        if self._pep_key:
            pep_raw = sii.get(self._pep_key)
            pep = float(pep_raw) if pep_raw is not None else None

        rank_value = sii.get("rank")
        rank = int(rank_value) if rank_value is not None else None

        psm = PSM(
            peptidoform=peptidoform,
            spectrum_id=psm_spectrum_id,
            run=run,
            is_decoy=is_decoy,
            score=score,
            qvalue=qvalue,
            pep=pep,
            precursor_mz=precursor_mz,
            retention_time=rt,
            ion_mobility=ion_mobility,
            protein_list=protein_list,
            rank=rank,
            source=self._source,
            provenance_data={"mzid_filename": str(self.filename)},
            metadata=metadata,
        )
        return psm

    def _get_non_metadata_keys(self, keys: list[str]) -> None:
        """Gather all the keys at PSM-level that should not be written to metadata."""
        # All keys required to create PSM object
        default_keys = [
            "chargeState",
            "rank",
            "PeptideSequence",
            "experimentalMassToCharge",
            "PeptideEvidenceRef",
            "Modification",
        ]
        # Get the score key and add to default keys
        if not self._score_key:
            self._score_key = self._infer_score_name(keys)
        if self._score_key:
            default_keys.append(self._score_key)
        else:
            logger.warning(
                "No known score metric found in mzIdentML file. Scores will be set to None."
            )

        # Get the q-value key and add to default keys
        self._qvalue_key = self._infer_qvalue_name(keys)
        if self._qvalue_key:
            default_keys.append(self._qvalue_key)

        # Get the PEP key and add to default keys
        self._pep_key = self._infer_pep_name(keys)
        if self._pep_key:
            default_keys.append(self._pep_key)

        # Get retention time key
        for rt_key in ["retention time", "scan start time"]:
            if rt_key in keys:
                self._rt_key = rt_key
                default_keys.append(rt_key)
                break

        # Keys that are not necessary for metadata
        self._non_metadata_keys.extend(default_keys)

    def _get_toplevel_non_metadata_keys(self, keys: list[str]) -> None:
        """Identify spectrum-level keys that should not be written to PSM metadata."""
        # Check if RT is encoded in spectrum metadata
        for key in ["retention time", "scan start time"]:
            if key in keys:
                self._spectrum_rt_key = key
                self._non_metadata_keys.append(key)
                break

        # Check if ion mobility is encoded in spectrum metadata
        for im_key in ["inverse reduced ion mobility"]:
            if im_key in keys:
                self._im_key = im_key
                self._non_metadata_keys.append(im_key)
                break

    @staticmethod
    def _infer_score_name(keys: list[str]) -> str | None:
        """Infer search engine score name from available PSM keys."""
        lower_keys = {key.lower(): key for key in keys}
        for score in STANDARD_SEARCHENGINE_SCORES:
            if score in lower_keys:
                return lower_keys[score]
        return None

    @staticmethod
    def _infer_qvalue_name(keys: list[str]) -> str | None:
        """Infer q-value field name from available PSM keys."""
        for qvalue in Q_VALUE_TERMS:
            if qvalue in keys:
                return qvalue
        return None

    @staticmethod
    def _infer_pep_name(keys: list[str]) -> str | None:
        """Infer PEP (Posterior Error Probability) field name from available PSM keys."""
        for pep in PEP_TERMS:
            if pep in keys:
                return pep
        return None


class MzidWriter(WriterBase):
    """Writer for MzIdentMl files."""

    def __init__(
        self,
        filename: str | Path,
        *args: Any,
        show_progressbar: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Writer for mzIdentML PSM files.

        Parameters
        ----------
        filename: str, Pathlib.Path
            Path to PSM file.
        *args
            Additional positional argument passed to parent class.
        show_progressbar: bool, optional
            Show progress bar for conversion process. (default: False)
        **kwargs
            Additional keyword arguments passed to parent class.

        Notes
        -----
        - Unlike other psm_utils.io writer classes, :py:class:`MzidWriter` does not support writing
          a single PSM to a file with the :py:meth:`write_psm` method. Only writing a full PSMList
          to a file at once with the :py:meth:`write_file` method is currently supported.
        - While not required according to the mzIdentML specification document (v1.1.1), the
          retention time is written as cvParam ``retention time`` to the
          SpectrumIdentificationItem element. As the actual unit is not known in psm_utils,
          the unit is written as seconds.
        - As the actual PSM score type is not known in psm_utils, the score is written as cvParam
          ``MS:1001153`` to the SpectrumIdentificationItem element.

        """
        super().__init__(filename, *args, **kwargs)
        self.show_progressbar = show_progressbar
        self._writer = None

    def __enter__(self) -> MzidWriter:
        """Open file for writing and return self."""
        return self

    def __exit__(self, *args, **kwargs) -> None:
        """Close file and writer."""
        pass

    def write_psm(self, psm: PSM):
        """
        Write a single PSM to the PSM file.

        This method is currently not supported (see Notes).

        Raises
        ------
        NotImplementedError
            MzidWriter currently does not support write_psm.

        """
        raise NotImplementedError("MzidWriter currently does not support write_psm.")

    def write_file(self, psm_list: PSMList) -> None:
        """Write entire PSMList to mzIdentML file."""
        file = open(self.filename, "wb")
        with Progress(disable=not self.show_progressbar) as progress:
            with MzIdentMLWriter(file, close=True) as writer:
                writer.controlled_vocabularies()
                writer.provenance(
                    software={
                        "name": "psm_utils",
                        "uri": "https://github.com/compomics/psm_utils",
                        "version": __version__,
                    }
                )
                writer.register("SpectraData", 1)
                writer.register("SearchDatabase", 1)
                writer.register("SpectrumIdentificationList", 1)
                writer.register("SpectrumIdentificationProtocol", 1)

                proteins = set()
                peptide_ids = set()
                peptide_evidence_ids = set()

                proteins = {
                    prot
                    for prot_list in list(psm_list["protein_list"])
                    if prot_list
                    for prot in prot_list
                }

                spec_id_dict = psm_list.get_psm_dict()
                task1 = progress.add_task("[cyan]Writing Proteins to mzid", total=len(proteins))
                task2 = progress.add_task(
                    "[cyan]Writing Peptide and PeptideEvidence items",
                    total=len(psm_list),
                )
                task3 = progress.add_task(
                    "[cyan]Writing SpectrumIdentificationResults",
                    total=len(psm_list),
                )

                with writer.sequence_collection():
                    for prot in proteins:
                        writer.write_db_sequence(prot, None, id=prot, params=[])
                        progress.update(task1, advance=1)
                    for psm in psm_list:
                        peptide = psm["peptidoform"]
                        if peptide not in peptide_ids:
                            writer.write_peptide(**self._create_peptide_object(peptide))
                            peptide_ids.add(peptide)

                        if psm["protein_list"]:
                            for protein in psm["protein_list"]:
                                peptide_evidence_id = (
                                    f"PeptideEvidence_{peptide.proforma}_{protein}"
                                )
                                if peptide_evidence_id not in peptide_evidence_ids:
                                    peptide_evidence_ids.add(peptide_evidence_id)
                                    writer.write_peptide_evidence(
                                        peptide_id="Peptide_" + peptide.proforma,
                                        db_sequence_id=protein,
                                        id=peptide_evidence_id,
                                        start_position=None,
                                        end_position=None,
                                        is_decoy=psm["is_decoy"],
                                    )
                        progress.update(task2, advance=1)
                with writer.analysis_collection():
                    writer.SpectrumIdentification([1], [1]).write(writer)

                with writer.analysis_protocol_collection():
                    writer.spectrum_identification_protocol()  # build without?

                with writer.data_collection():
                    spectra_data, spectra_data_id_dict = self._transform_spectra_data(
                        spec_id_dict=spec_id_dict
                    )
                    writer.inputs(
                        source_files=[],
                        # # if fasta file is given, we can parse here and add protein information
                        # search_databases=transform_search_database(),
                        spectra_data=spectra_data,
                    )

                with writer.analysis_data():
                    with writer.spectrum_identification_list(id=1):
                        for collection in spec_id_dict.keys():
                            for run in spec_id_dict[collection].keys():
                                spectra_data_id = spectra_data_id_dict[
                                    "/".join(filter(None, [collection, run]))
                                ]
                                for spec_id in spec_id_dict[collection][run].keys():
                                    identified_psms = spec_id_dict[collection][run][spec_id]
                                    writer.write_spectrum_identification_result(
                                        **self._transform_spectrum_identification_result(
                                            spec_id, identified_psms, spectra_data_id
                                        )
                                    )
                                    progress.update(
                                        task3,
                                        advance=len(spec_id_dict[collection][run][spec_id]),
                                    )

    @staticmethod
    def _create_peptide_object(peptidoform: Peptidoform) -> dict[str, Any]:
        """Create mzid peptide object from Peptidoform."""

        def parse_modifications(
            modifications: list[Any] | None, location: int
        ) -> list[dict[str, Any]]:
            modification_list = []
            if modifications:
                for mod in modifications:
                    try:
                        modification_list.append(
                            {
                                "location": location,
                                "name": mod.name,  # type: ignore[attr-defined]
                                "monoisotopic_mass_delta": mod.mass,  # type: ignore[attr-defined]
                            }
                        )
                    except AttributeError:
                        modification_list.append(
                            {
                                "location": location,
                                "monoisotopic_mass_delta": mod.mass,  # type: ignore[attr-defined]
                            }
                        )
            return modification_list

        # Parse modifications
        modifications = []
        for loc, (aa, mods) in enumerate(peptidoform.parsed_sequence, start=1):
            modifications.extend(parse_modifications(mods, loc))
        modifications.extend(parse_modifications(peptidoform.properties.get("n_term"), 0))
        modifications.extend(
            parse_modifications(
                peptidoform.properties.get("c_term"), len(peptidoform.sequence) + 1
            )
        )

        peptide_object = {
            "id": "Peptide_" + peptidoform.proforma,
            "peptide_sequence": peptidoform.sequence,
            "modifications": modifications,
        }

        return peptide_object

    def _transform_search_database(self) -> dict[str, Any]:
        """Create mzid database object."""
        # TODO: Create this and link with protein object when fasta file is provided
        return {
            "file_format": "fasta format",
            "name": "",
            "id": 1,
            "location": "",
            "params": [],
        }

    @staticmethod
    def _transform_spectra_data(
        spec_id_dict: dict[str, Any],
    ) -> tuple[list[dict[str, Any]], dict[str, int]]:
        """Get all unique spectra data from PSMList spectrum id dict."""
        collection_run_id_dict: dict[str, int] = {}
        spectra_data = []
        i = 1
        for collection in spec_id_dict:
            for run in spec_id_dict[collection]:
                collection_run_id = "/".join(filter(None, [collection, run]))
                if collection_run_id not in collection_run_id_dict:
                    collection_run_id_dict[collection_run_id] = i
                    spectra_data_object = {
                        "id": i,
                        "location": collection_run_id,
                        "spectrum_id_format": "multiple peak list nativeID format",
                        # 'file_format': #TODO can we infer this?
                    }
                    spectra_data.append(spectra_data_object)
                    i += 1
        return spectra_data, collection_run_id_dict

    @staticmethod
    def _transform_spectrum_identification_item(
        candidate_psm: dict[str, Any],
    ) -> list[dict[str, Any]]:
        """Create SpectrumIdentificationItem for each candidate PSM."""
        peptide = candidate_psm["peptidoform"].proforma
        params = list(candidate_psm["metadata"].items()) if candidate_psm.get("metadata") else []
        params = [{k: v} for k, v in params]

        for key, label, unit in [
            ("retention_time", "retention time", "second"),
            ("qvalue", "PSM-level q-value", None),
            ("pep", "PSM-level local FDR", None),
        ]:
            value = candidate_psm.get(key)
            if value is not None:
                param = {"name": label, "value": value}
                if unit:
                    param["unit_name"] = unit
                params.append(param)

        candidate_psm_dict = {
            "charge_state": candidate_psm["peptidoform"].precursor_charge,
            "peptide_id": f"Peptide_{peptide}",
            "score": {"search engine specific score": candidate_psm["score"]},
            "experimental_mass_to_charge": candidate_psm["precursor_mz"],
            "calculated_mass_to_charge": candidate_psm["peptidoform"].theoretical_mz,
            "rank": candidate_psm["rank"],
            "params": params,
        }
        items = []
        if candidate_psm["protein_list"]:
            for prot_acc in candidate_psm["protein_list"]:
                protein_specific_items = {
                    "peptide_evidence_id": f"PeptideEvidence_{peptide}_{prot_acc}",
                    "id": f"SII_{candidate_psm['spectrum_id']}_{peptide}_{prot_acc}",
                }
                items.append(dict(candidate_psm_dict, **protein_specific_items))
        return items

    def _transform_spectrum_identification_result(
        self, spec_id: str, identified_psms: list[dict[str, Any]], spectra_data_id: int
    ) -> dict[str, Any]:
        """Create mzid SpectrumIdentificationResult object for spectrum PSMs."""
        spectrum_id_result = {
            "id": f"SIR_{spec_id}",
            "spectrum_id": spec_id,
            "spectra_data_id": spectra_data_id,
        }
        identifications = []
        for candidate in identified_psms:
            identifications.extend(self._transform_spectrum_identification_item(candidate))
        spectrum_id_result["identifications"] = identifications
        return spectrum_id_result


class UnknownMzidScore(PSMUtilsIOException):
    """No known score metric found in mzIdentML file."""
