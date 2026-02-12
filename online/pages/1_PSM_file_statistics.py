"""psm_utils Streamlit-based web server."""

import gzip
import re
from tempfile import NamedTemporaryFile

import numpy as np
import streamlit as st
from _base import StreamlitPage
from _utils import fdr_plot, pp_plot, score_histogram

from psm_utils.io import READERS, _infer_filetype, read_file


class StreamlitPageStats(StreamlitPage):
    """Statistics page for psm_utils online Streamlit web server."""

    def _main_page(self):
        st.markdown(
            """
            ## PSM file statistics
            """
        )

        if self._input_form():
            try:
                self._read_file()
                self._prepare_psms()
            except Exception as e:
                st.error("Error occurred while reading PSM file.")
                st.exception(e)
            else:
                self._show_results()

    def _input_form(self):
        with st.form(key="main_form"):
            self.state["input_file"] = st.file_uploader(
                label="Input PSM file",
                accept_multiple_files=False,
                help=(
                    "PSM file. See https://psm_utils.readthedocs.io for all supported "
                    "file types. Uploaded files are limited to 200MB, although single "
                    "GZipped files are also supported (e.g., `msms.txt.gz`). See "
                    "https://github.com/compomics/psm_utils/tree/main/example_files "
                    "for examples."
                ),
            )
            self.state["input_filetype"] = st.selectbox(
                label="Input PSM file type",
                options=["infer"] + list(READERS.keys()),
                index=0,
            )
            self.state["use_example_file"] = st.checkbox(
                label="Use example PSM file",
            )

            with st.expander("Advanced options"):
                row = st.columns(3)
                self.state["decoy_pattern"] = row[0].text_input(
                    "Decoy pattern",
                    help=(
                        """
                        Decoy protein entries are commonly marked with a prefix or
                        suffix, e.g. `DECOY_`, or `_REVERSED`. If the target/decoy state
                        of a PSM is not explicitly encoded in the PSM file, this setting
                        can be used to find decoy PSMs by trying to match a regular
                        expression pattern to the protein names. Some patterns that
                        might work are `^DECOY_` or `_REVERSED$`.
                        """
                    ),
                )
                self.state["fdr_threshold"] = row[1].number_input(
                    label="FDR threshold",
                    min_value=0.0,
                    max_value=1.0,
                    step=0.001,
                    value=0.01,
                    format="%f",
                )
                self.state["percolator_score_column"] = row[2].text_input(
                    "Percolator Tab: Score column",
                    help=(
                        """
                        In Percolator Tab PIN files, the name of the score column is not
                        predefined. Provide the correct column name to extract PSM
                        scores from a PIN file.
                        """
                    ),
                )
                row = st.columns(2)
                self.state["reverse"] = row[0].radio(
                    "Score type: order",
                    options=[True, False],
                    format_func=lambda x: (
                        "Higher score is better" if x else "Lower score is better"
                    ),
                )
                self.state["log_scale"] = row[1].radio(
                    "Score type: scale",
                    options=[False, True],
                    format_func=lambda x: (
                        "Logarithmic scale (e.g., e-value)"
                        if x
                        else "Linear scale (e.g., Andromeda score)"
                    ),
                    help=(
                        """
                        Some search engine scores, mostly e-value-like scores, require a
                        logarithmic transformation for plotting. Usually, these scores
                        also require the "Lower score is better" option.
                        """
                    ),
                )

            submitted = st.form_submit_button("Upload")

        if submitted:
            if not (self.state["input_file"] or self.state["use_example_file"]):
                st.error(
                    "Input PSM file required. Upload a file or select 'Use example PSM file'."
                )
            else:
                return True

    def _read_file(self):
        with st.spinner("Reading PSM file..."):
            # Reading file
            if self.state["use_example_file"]:
                psm_list = read_file(
                    "online/example_data/LFQ_Orbitrap_DDA_Ecoli_01_all.tsv",
                    filetype="tsv",
                )
            else:
                # Infer filetype if required
                if st.session_state["input_filetype"] == "infer":
                    st.session_state["input_filetype"] = _infer_filetype(
                        re.sub(
                            r".gz$",
                            "",
                            st.session_state["input_file"].name,
                            flags=re.IGNORECASE,
                        )
                    )
                # Write file to disk for psm_utils; then read
                with NamedTemporaryFile(mode="wb", delete=False) as tmp_file:
                    if self.state["input_file"].name.lower().endswith(".gz"):
                        tmp_file.write(gzip.decompress(self.state["input_file"].getvalue()))
                    else:
                        tmp_file.write(self.state["input_file"].getvalue())
                    tmp_file.flush()
                    tmp_file.close()
                    psm_list = read_file(
                        filename=tmp_file.name,
                        filetype=self.state["input_filetype"],
                        score_column=self.state["percolator_score_column"],
                    )
        if self.state["decoy_pattern"]:
            psm_list.find_decoys(self.state["decoy_pattern"])

        self.state["psm_list"] = psm_list

    def _prepare_psms(self):
        """Calculate q-values if needed; generate psm_df."""
        psm_list = self.state["psm_list"]

        # Transform logarithmic score, if needed
        if self.state["log_scale"]:
            psm_list["score"] = np.log2(psm_list["score"])

        # Warn if some but only few decoys are present
        percent_decoys = np.count_nonzero(psm_list["is_decoy"]) / len(psm_list)
        if 0 < percent_decoys < 0.1:
            st.warning(
                f"""
                An unusually low percentage of decoy PSMs was found
                ({percent_decoys:.2%}). Are you sure that the PSM file was not already
                FDR-filtered?
                """
            )

        # If no q-values, try to calculate
        if (psm_list["qvalue"] == None).any():  # noqa: E711
            # If no decoys, display error
            if percent_decoys == 0.0:
                st.error(
                    """
                    No decoys were found in the PSM file and not all PSMs have q-values
                    assigned. The FDR could therefore not be calculated. Are you sure
                    that both target and decoy PSMs were returned by the search engine?
                    If there should be decoys present in the PSM file, you might want
                    to set the advanced 'Decoy pattern' option.
                    """
                )
                file_state = ["no_qvalues", "no_decoys"]
            else:
                st.warning(
                    """
                    Not all PSMs have q-values assigned, so q-values will be calculated
                    based on target and decoy PSM scores. Please ensure that all decoy
                    PSMs are present in the PSM file to prevent an incorrect
                    FDR estimation.
                    """
                )
                try:
                    psm_list.calculate_qvalues(self.state["reverse"])
                except Exception as e:
                    st.exception(e)
                    file_state = ["no_qvalues", "decoys"]
                else:
                    file_state = ["qvalues", "decoys"]
        else:
            if percent_decoys == 0.0:
                file_state = ["qvalues", "no_decoys"]
            else:
                file_state = ["qvalues", "decoys"]

        self.state["file_state"] = file_state
        self.state["psm_df"] = psm_list.to_dataframe()

    def _show_results(self):
        # General stats
        st.markdown(
            """
            ## Results
            ### Overall statistics
            Total number of items in the PSM file.
            """
        )
        psm_list = self.state["psm_list"]
        psm_df = self.state["psm_df"]

        n_collections = psm_df["collection"].unique().shape[0]
        n_runs = psm_df[["run", "collection"]].drop_duplicates().shape[0]
        n_spectra = psm_df[["spectrum_id", "run", "collection"]].drop_duplicates().shape[0]
        n_psms = psm_df.shape[0]
        n_peptidoforms = psm_df["peptidoform"].apply(lambda x: x.proforma).unique().shape[0]
        percent_decoys = np.count_nonzero(psm_list["is_decoy"]) / len(psm_list)

        row_1 = st.columns(3)
        row_1[0].metric(label="Collections", value=n_collections)
        row_1[1].metric(label="Runs", value=n_runs)
        row_1[2].metric(label="Spectra", value=n_spectra)
        row_2 = st.columns(3)
        row_2[0].metric(label="PSMs", value=n_psms)
        row_2[1].metric(label="Peptidoforms", value=n_peptidoforms)
        row_2[2].metric(label="Decoy PSMs", value=f"{percent_decoys:.2%}")

        # FDR-filtered stats
        st.markdown(
            """
            ### FDR-filtered statistics
            Number of identifications filtered at the FDR threshold as selected in
            the input form.
            """
        )
        if "no_qvalues" in self.state["file_state"]:
            st.error("PSM q-values are required to filter on FDR.")
        else:
            psm_df_filtered = psm_df[psm_df["qvalue"] <= self.state["fdr_threshold"]]
            n_spectra = (
                psm_df_filtered[["spectrum_id", "run", "collection"]].drop_duplicates().shape[0]
            )
            n_psms = psm_df_filtered.shape[0]
            n_peptides = psm_df["peptidoform"].apply(lambda x: x.sequence).unique().shape[0]
            n_peptidoforms = (
                psm_df_filtered["peptidoform"].apply(lambda x: x.proforma).unique().shape[0]
            )

            row_3 = st.columns(4)
            row_3[0].metric(label="Spectra", value=n_spectra)
            row_3[1].metric(label="PSMs", value=n_psms)
            row_3[2].metric(label="Peptides", value=n_peptides)
            row_3[3].metric(label="Peptidoforms", value=n_peptidoforms)

        # Plots
        st.markdown(
            """
            ### Target-decoy diagnostic plots
            #### Score histogram
            The score histogram shows the score distribution for both target and
            decoy PSMs.
            """
        )
        try:
            st.plotly_chart(score_histogram(psm_df), use_container_width=True)
        except Exception as e:
            st.error(e)

        st.markdown(
            """
            #### Percentile-percentile plot
            The percentile-percentile (PP) plot shows the empirical cumulative
            distribution function (ECDF) of the target distribution in function of
            the ECDF of the decoy distribution. In the context of peptide
            identification, it can be used to assess the quality of decoy PSMs and
            their capacity to help in correctly estimating the false discovery rate.

            Ideally, the PP-plot should follow a straight diagonal line up until the
            end of the decoy distribution (right-hand side of the plot), where the
            line turns vertically upwards. This means that the decoy distribution
            perfectly aligns with the first part of the target distribution (the
            low-scoring and presumably bad target PSMs) and therefore correctly
            models the bad target PSMs. This diagonal line matches the ratio of
            the number of decoy to the number of target PSMs ($\widehat{\pi}_0$).

            More information on this type of diagnostic plot can be found at
            [statomics.github.io/TargetDecoy](https://statomics.github.io/TargetDecoy/articles/TargetDecoy.html).

            """
        )
        if "no_decoys" in self.state["file_state"]:
            st.error("Decoy PSMs are required to generate this plot.")
        else:
            try:
                st.plotly_chart(pp_plot(psm_df), use_container_width=True)
            except Exception as e:
                st.exception(e)

        st.markdown(
            """
            #### False discovery rate plot
            This plot shows the number of identified target PSMs in function of the
            FDR threshold. The plot starts at the top-right corner with
            the total number of PSMs in the dataset (no FDR filtering). As the FDR
            threshold becomes more stringent (towards the left of the x-axis), the
            number of identified target PSMs goes down. The red vertical line
            indicates the FDR threshold as selected in the input form.
            """
        )
        if "no_qvalues" in self.state["file_state"]:
            st.error("PSM q-values are required to generate this plot.")
        else:
            try:
                st.plotly_chart(
                    fdr_plot(psm_df, self.state["fdr_threshold"]),
                    use_container_width=True,
                )
            except Exception as e:
                st.exception(e)


if __name__ == "__main__":
    StreamlitPageStats()
