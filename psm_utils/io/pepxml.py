"""Interface with TPP pepXML PSM files."""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from pyteomics import mass, pepxml, proforma  # type: ignore[import]

from psm_utils.io._base_classes import ReaderBase
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.utils import mass_to_mz

logger = logging.getLogger(__name__)

STANDARD_SEARCHENGINE_SCORES = [
    "expect",
    "EValue",
    "Evalue",
    "SpecEValue",
    "xcorr",  # Fallback if no e-value is present
    "delta_dot",  # SpectraST
    "mzFidelity",
]


class PepXMLReader(ReaderBase):
    """Reader for pepXML PSM files."""

    def __init__(
        self, filename: str | Path, *args: Any, score_key: str | None = None, **kwargs: Any
    ) -> None:
        """
        Reader for pepXML PSM files.

        Parameters
        ----------
        filename
            Path to PSM file.
        *args
            Additional positional arguments passed to parent class.
        score_key
            Name of the score metric to use as PSM score. If not provided, the score metric is
            inferred from a list of known search engine scores.
        **kwargs
            Additional keyword arguments passed to parent class.

        """
        super().__init__(filename, *args, **kwargs)
        self.score_key = score_key or self._infer_score_name()

    def __iter__(self) -> Iterator[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with pepxml.read(str(self.filename)) as reader:
            for spectrum_query in reader:
                if "search_hit" not in spectrum_query:
                    continue
                for search_hit in spectrum_query["search_hit"]:
                    yield self._parse_psm(spectrum_query, search_hit)

    def _infer_score_name(self) -> str | None:
        """Infer the score from the list of known PSM scores."""
        # Get scores from first PSM
        with pepxml.read(str(self.filename)) as reader:
            for spectrum_query in reader:
                score_keys = spectrum_query["search_hit"][0]["search_score"].keys()
                break
            else:
                score_keys = []

        # Infer score name
        if not score_keys:
            logger.warning("No pepXML scores found.")
            return None

        for score in STANDARD_SEARCHENGINE_SCORES:  # Check for known scores
            if score in score_keys:
                logger.debug(f"Using known pepXML score `{score}`.")
                return score

        # Default to the first one if nothing found
        logger.warning(f"No known pepXML scores found. Defaulting to `{score_keys[0]}`.")
        return score_keys[0]

    @staticmethod
    def _parse_peptidoform(
        peptide: str, modifications: list[dict[str, Any]], charge: int | None = None
    ) -> Peptidoform:
        """Parse pepXML peptide to Peptidoform."""
        modifications_dict: dict[int, list[Any]] = defaultdict(list)
        n_term: list[Any] = []
        c_term: list[Any] = []

        for mod in modifications:
            # Round mass modification to 6 decimal places, precision from UniMod
            if mod["position"] == 0:
                mod_tag = proforma.process_tag_tokens(f"{mod['mass']:+.6f}")
                n_term.append(mod_tag)
            elif mod["position"] == len(peptide) + 1:
                mod_tag = proforma.process_tag_tokens(f"{mod['mass']:+.6f}")
                c_term.append(mod_tag)
            else:
                # Convert 1-based to 0-based position
                position = mod["position"] - 1
                # Sequence modifications are written as residue mass + modification mass
                mod_mass = mod["mass"] - mass.std_aa_mass[peptide[position]]
                mod_tag = proforma.process_tag_tokens(f"{mod_mass:+.6f}")
                modifications_dict[position].append(mod_tag)

        sequence = [(aa, modifications_dict[i] or None) for i, aa in enumerate(peptide)]
        properties = {
            "n_term": n_term,
            "c_term": c_term,
            "charge_state": proforma.ChargeState(charge) if charge else None,
            "unlocalized_modifications": [],
            "labile_modifications": [],
            "fixed_modifications": [],
            "intervals": [],
            "isotopes": [],
            "group_ids": [],
        }
        return Peptidoform(proforma.ProForma(sequence, properties))

    def _parse_psm(self, spectrum_query: dict[str, Any], search_hit: dict[str, Any]) -> PSM:
        """Parse pepXML PSM to PSM."""
        metadata = {
            "num_matched_ions": str(search_hit["num_matched_ions"]),
            "tot_num_ions": str(search_hit["tot_num_ions"]),
            "num_missed_cleavages": str(search_hit["num_missed_cleavages"]),
        }
        metadata.update(
            {
                f"search_score_{key.lower()}": str(search_hit["search_score"][key])
                for key in search_hit["search_score"]
            }
        )

        return PSM(
            peptidoform=self._parse_peptidoform(
                search_hit["peptide"],
                search_hit["modifications"],
                spectrum_query["assumed_charge"],
            ),
            spectrum_id=spectrum_query["spectrumNativeID"]
            if "spectrumNativeID" in spectrum_query
            else spectrum_query["spectrum"],
            run=None,
            collection=None,
            spectrum=None,
            is_decoy=None,
            score=search_hit["search_score"][self.score_key],
            qvalue=None,
            pep=None,
            precursor_mz=mass_to_mz(
                spectrum_query["precursor_neutral_mass"], spectrum_query["assumed_charge"]
            ),
            retention_time=spectrum_query.get("retention_time_sec"),
            ion_mobility=spectrum_query.get("ion_mobility"),
            protein_list=[p["protein"] for p in search_hit["proteins"]],
            rank=search_hit["hit_rank"],
            source=None,
            provenance_data={
                "pepxml_index": str(spectrum_query["index"]),
                "start_scan": str(spectrum_query["start_scan"]),
                "end_scan": str(spectrum_query["end_scan"]),
            },
            metadata=metadata,
            rescoring_features={},
        )
