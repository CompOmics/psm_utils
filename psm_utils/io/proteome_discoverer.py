"""
Reader for Proteome Discoverer MSF PSM files.

This module provides functionality to read PSM data from Proteome Discoverer MSF SQLite database
files.

The reader supports both target and decoy peptides, handles various modification types (amino acid
and terminal modifications), and extracts complete scoring information from the MSF database
structure.

Examples
--------
>>> from psm_utils.io.proteome_discoverer import MSFReader
>>> reader = MSFReader("results.msf")
>>> psm_list = reader.read_file()
>>> for psm in reader:
...     print(f"{psm.peptidoform} - Score: {psm.score}")

Notes
-----
MSF file versions 79, 53, and 8 are currently supported.

"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pyteomics.proforma as proforma  # type: ignore[import-untyped]
from sqlalchemy import create_engine, func, select
from sqlalchemy.orm import Session

import psm_utils.io._pd_msf_tables as msf
from psm_utils import PSM, Peptidoform
from psm_utils.io._base_classes import ReaderBase

logger = logging.getLogger(__name__)

COMPATIBLE_VERSIONS = [79, 53, 8]


class MSFReader(ReaderBase):
    """
    Reader for Proteome Discoverer MSF files.

    This reader processes SQLite-based MSF database files from Proteome Discoverer, extracting
    peptide-spectrum matches with complete modification information, scoring data, and metadata.
    Supports both target and decoy peptides.

    Examples
    --------
    >>> reader = MSFReader("experiment.msf")
    >>> psm_list = reader.read_file()
    >>> len(reader)  # Get total number of PSMs
    1234
    >>> for psm in reader:  # Iterate over all PSMs
    ...     if psm.qvalue and psm.qvalue < 0.01:
    ...         print(f"High-confidence PSM: {psm.peptidoform}")

    """

    def __init__(
        self,
        filename: str | Path,
        *args,
        **kwargs,
    ) -> None:
        """
        Initialize MSF reader with database connection and version validation.

        Parameters
        ----------
        filename
            Path to Proteome Discoverer MSF file.
        *args
            Additional positional arguments passed to parent class.
        **kwargs
            Additional keyword arguments passed to parent class.

        """
        super().__init__(filename, *args, **kwargs)

        self._engine = create_engine(f"sqlite:///{self.filename.as_posix()}")
        self._session = Session(self._engine)

        self._check_version()

    def __len__(self) -> int:
        """Return total number of PSMs in file."""
        peptide_count = (
            self._session.execute(select(func.count()).select_from(msf.Peptide)).scalar() or 0
        )
        decoy_count = (
            self._session.execute(select(func.count()).select_from(msf.PeptideDecoy)).scalar() or 0
        )
        return peptide_count + decoy_count

    def __iter__(self) -> Iterator[PSM]:
        """
        Iterate over file and return PSMs one-by-one.

        Yields
        ------
        PSM
            Individual PSM objects with complete modification and scoring information.

        """
        for is_decoy in [False, True]:
            modifications = self._get_modifications(is_decoy)
            terminal_modifications = self._get_terminal_modifications(is_decoy)
            protein_entries = self._get_protein_entries(is_decoy)
            main_score = self._get_main_score(is_decoy)
            secondary_scores = self._get_secondary_scores(is_decoy)

            for entry in self._iter_peptides(is_decoy):
                peptide = entry[0]  # First element is Peptide or PeptideDecoy
                peptide_id = peptide.PeptideID
                yield self._parse_entry(
                    entry,
                    modifications[peptide_id],
                    terminal_modifications[peptide_id],
                    protein_entries[peptide_id],
                    main_score[peptide_id],
                    secondary_scores[peptide_id],
                    is_decoy,
                )

    def __enter__(self) -> ReaderBase:
        """Enter context manager for MSFReader."""
        return super().__enter__()

    def __exit__(self, *args, **kwargs) -> None:
        """Exit context manager for MSFReader."""
        self._session.close()
        return super().__exit__(*args, **kwargs)

    def _check_version(self) -> None:
        """Check MSF file version compatibility."""
        first_result = self._session.execute(select(msf.SchemaInfo.Version)).first()
        if first_result is None:
            logger.warning(
                "MSF file does not contain version information. "
                "Assuming it is compatible with this reader."
            )
            return None
        version = first_result[0]
        if version not in COMPATIBLE_VERSIONS:
            logger.warning(
                f"MSF file version {version} might not be compatible with this reader. "
                f"Checked versions are: {COMPATIBLE_VERSIONS}."
            )

    def _iter_peptides(self, is_decoy: bool) -> Iterator[Any]:
        """
        Iterate over peptides in MSF file.

        Parameters
        ----------
        is_decoy : bool
            Whether to iterate over decoy peptides instead of target peptides.

        Yields
        ------
        Any
            SQLAlchemy Row object containing joined Peptide, SpectrumHeader, MassPeak, and
            FileInfo data. The Row object has attributes like row[0] (Peptide/PeptideDecoy),
            row[1] (SpectrumHeader), row[2] (MassPeak), and row[3] (FileInfo).

        Notes
        -----
        This method performs a complex join across multiple MSF database tables to gather
        all necessary information for PSM construction. The returned rows contain all
        spectral and identification metadata needed for downstream processing.

        """
        # Select appropriate peptide table based on decoy flag
        peptide_table = msf.PeptideDecoy if is_decoy else msf.Peptide

        # Build and execute query - same structure for both target and decoy
        stmt = (
            select(peptide_table, msf.SpectrumHeader, msf.MassPeak, msf.FileInfo)
            .select_from(peptide_table)
            .join(msf.SpectrumHeader, peptide_table.SpectrumID == msf.SpectrumHeader.SpectrumID)
            .join(msf.MassPeak, msf.MassPeak.MassPeakID == msf.SpectrumHeader.MassPeakID)
            .join(msf.FileInfo, msf.FileInfo.FileID == msf.MassPeak.FileID)
        )

        yield from self._session.execute(stmt)

    def _get_modifications(self, is_decoy: bool) -> dict[int, list[tuple[int, int]]]:
        """Get amino acid modifications per peptide ID."""
        PeptidesAminoAcidModification = (
            msf.PeptidesAminoAcidModificationsDecoy
            if is_decoy
            else msf.PeptidesAminoAcidModification
        )
        stmt = (
            select(
                PeptidesAminoAcidModification.PeptideID,
                PeptidesAminoAcidModification.Position,
                msf.AminoAcidModification.UnimodAccession,
            )
            .select_from(PeptidesAminoAcidModification)
            .join(
                msf.AminoAcidModification,
                PeptidesAminoAcidModification.AminoAcidModificationID
                == msf.AminoAcidModification.AminoAcidModificationID,
            )
        )
        modifications_by_peptide: dict[int, list[tuple[int, int]]] = defaultdict(list)
        for peptide_id, position, unimod_accession in self._session.execute(stmt):
            modifications_by_peptide[peptide_id].append((position, unimod_accession))

        return modifications_by_peptide

    def _get_terminal_modifications(self, is_decoy: bool) -> dict[int, list[tuple[int, int]]]:
        """Get terminal modifications per peptide ID."""
        PeptidesTerminalModification = (
            msf.PeptidesTerminalModification if is_decoy else msf.PeptidesTerminalModificationDecoy
        )
        stmt = (
            select(
                PeptidesTerminalModification.PeptideID,
                msf.AminoAcidModification.PositionType,
                msf.AminoAcidModification.UnimodAccession,
            )
            .select_from(msf.AminoAcidModification)
            .join(
                PeptidesTerminalModification,
                PeptidesTerminalModification.TerminalModificationID
                == msf.AminoAcidModification.AminoAcidModificationID,
            )
        )
        terminal_modifications: dict[int, list[tuple[int, int]]] = defaultdict(list)
        for peptide_id, position_type, unimod_accession in self._session.execute(stmt):
            terminal_modifications[peptide_id].append((position_type, unimod_accession))
        return terminal_modifications

    def _get_protein_entries(self, is_decoy: bool) -> dict[int, list[str]]:
        """Get protein descriptions per peptide ID."""
        PeptidesProtein = msf.PeptidesProteinDecoy if is_decoy else msf.PeptidesProtein
        stmt = (
            select(PeptidesProtein.PeptideID, msf.ProteinAnnotation.Description)
            .select_from(PeptidesProtein)
            .join(
                msf.ProteinAnnotation,
                PeptidesProtein.ProteinID == msf.ProteinAnnotation.ProteinID,
            )
        )
        proteins: dict[int, list[str]] = defaultdict(list)
        for peptide_id, description in self._session.execute(stmt):
            proteins[peptide_id].append(re.sub(r"^>", "", description))
        return proteins

    def _get_main_score(self, is_decoy: bool) -> dict[int, tuple[float, str]]:
        """Get main score and name per peptide ID."""
        PeptideScore = msf.PeptideScoreDecoy if is_decoy else msf.PeptideScore
        stmt = (
            select(
                PeptideScore.PeptideID, PeptideScore.ScoreValue, msf.ProcessingNodeScore.ScoreName
            )
            .select_from(PeptideScore)
            .join(
                msf.ProcessingNodeScore,
                msf.ProcessingNodeScore.ScoreID == PeptideScore.ScoreID,
            )
            .filter(msf.ProcessingNodeScore.IsMainScore == True)  # noqa: E712
        )
        scores: dict[int, tuple[float, str]] = {}
        for peptide_id, score_value, score_name in self._session.execute(stmt):
            scores[peptide_id] = (score_value, score_name)
        return scores

    def _get_secondary_scores(self, is_decoy: bool) -> dict[int, dict[str, float]]:
        """Get secondary scores per peptide ID."""
        PeptideScore = msf.PeptideScoreDecoy if is_decoy else msf.PeptideScore
        stmt = (
            select(
                PeptideScore.PeptideID, PeptideScore.ScoreValue, msf.ProcessingNodeScore.ScoreName
            )
            .select_from(PeptideScore)
            .join(
                msf.ProcessingNodeScore,
                msf.ProcessingNodeScore.ScoreID == PeptideScore.ScoreID,
            )
            .filter(msf.ProcessingNodeScore.IsMainScore == False)  # noqa: E712
        )
        scores: dict[int, dict[str, float]] = defaultdict(dict)
        for peptide_id, score_value, score_name in self._session.execute(stmt):
            scores[peptide_id][score_name] = score_value
        return scores

    def _compile_peptidoform(
        self,
        sequence: str,
        charge: int,
        modifications: list[tuple[int, int]],
        terminal_modifications: list[tuple[int, int]],
    ) -> Peptidoform:
        """
        Compile a peptidoform from a sequence, charge, and list of (terminal) modifications.

        Parameters
        ----------
        sequence
            The stripped sequence of the peptidoform.
        charge
            Precursor charge.
        modifications
            List of tuples of the form (position, unimod identifier).
        terminal_modifications
            List of tuples of the form (position type, unimod identifier).

        Notes
        -----
        The position type is either 1 (Any N-term), 2 (Any C-term), 3 (Protein N-term), or 4
        (Protein C-term). Position type 0 (Anywhere) should not be present in the
        terminal_modifications list.

        """
        modifications_dict = defaultdict(list)
        for position, unimod_id in modifications:
            modifications_dict[position].append(proforma.process_tag_tokens(f"U:{unimod_id}"))

        n_term = [
            proforma.process_tag_tokens(f"U:{unimod_id}")
            for position_type, unimod_id in terminal_modifications
            if position_type in [1, 3]  # Position types 'Any N-term' or 'Protein N-term'
        ]
        c_term = [
            proforma.process_tag_tokens(f"U:{unimod_id}")
            for position_type, unimod_id in terminal_modifications
            if position_type in [2, 4]  # Position types 'Any C-term' or 'Protein C-term'
        ]

        parsed_sequence = [(aa, modifications_dict[i] or None) for i, aa in enumerate(sequence)]
        properties = {
            "n_term": n_term,
            "c_term": c_term,
            "charge_state": proforma.ChargeState(charge),
            "unlocalized_modifications": [],
            "labile_modifications": [],
            "fixed_modifications": [],
            "intervals": [],
            "isotopes": [],
            "group_ids": [],
        }

        return Peptidoform(proforma.ProForma(parsed_sequence, properties))

    def _parse_entry(
        self,
        entry: Any,  # SQLAlchemy Row[tuple[Peptide|PeptideDecoy, SpectrumHeader, MassPeak, FileInfo]]
        modifications: list[tuple[int, int]],
        terminal_modifications: list[tuple[int, int]],
        protein_entries: list[str],
        main_score: tuple[float, str],
        secondary_scores: dict[str, float],
        is_decoy: bool,
    ) -> PSM:
        """
        Parse an entry from the MSF file into a PSM object.

        Parameters
        ----------
        entry : Any
            SQLAlchemy Row object containing joined peptide, spectrum, and file information.
            Accessed by index: entry[0] (Peptide/PeptideDecoy), entry[1] (SpectrumHeader),
            entry[2] (MassPeak), entry[3] (FileInfo).
        modifications : list[tuple[int, int]]
            List of tuples containing (position, UNIMOD accession) for amino acid modifications.
        terminal_modifications : list[tuple[int, int]]
            List of tuples containing (position_type, UNIMOD accession) for terminal modifications.
        protein_entries : list[str]
            List of protein descriptions associated with this peptide.
        main_score : tuple[float, str]
            Tuple containing (score_value, score_name) for the main search engine score.
        secondary_scores : dict[str, float]
            Dictionary mapping score names to values for secondary scores.
        is_decoy : bool
            Whether this PSM is from a decoy search.

        Returns
        -------
        PSM
            Complete PSM object with all available metadata and scoring information.

        Notes
        -----
        This method constructs a complete PSM object by:
        - Creating a peptidoform from sequence and modifications
        - Extracting spectrum identification and precursor information
        - Including all available scoring metrics
        - Adding proteome discoverer-specific metadata

        """
        peptide = entry[0]  # First element is Peptide or PeptideDecoy
        spectrum_header = entry[1]  # Second element is SpectrumHeader
        mass_peak = entry[2]  # Third element is MassPeak
        file_info = entry[3]  # Fourth element is FileInfo

        return PSM(
            peptidoform=self._compile_peptidoform(
                peptide.Sequence,
                spectrum_header.Charge,
                modifications,
                terminal_modifications,
            ),
            spectrum_id=spectrum_header.LastScan,
            run=Path(file_info.FileName).stem,
            is_decoy=is_decoy,
            score=main_score[0],
            qvalue=None,
            pep=None,
            precursor_mz=mass_peak.Mass,
            retention_time=spectrum_header.RetentionTime,
            ion_mobility=None,
            protein_list=protein_entries,
            rank=peptide.SearchEngineRank,
            source="proteome_discoverer",
            provenance_data={
                "scan_numbers": spectrum_header.ScanNumbers,
            },
            metadata={
                "ms1_intensity": str(mass_peak.Intensity),
                "ms1_percent_isolation_interference": str(mass_peak.PercentIsolationInterference),
                "ms1_ion_inject_time": str(mass_peak.IonInjectTime),
                "main_score_name": main_score[1],
                **secondary_scores,
            },
            rescoring_features={
                "missed_cleavages": peptide.MissedCleavages,
                "total_ions_count": peptide.TotalIonsCount,
                "matched_ions_count": peptide.MatchedIonsCount,
            },
        )
