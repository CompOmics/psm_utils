# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.5.0] - 2025-10-27

### Added

- ✨ `io`: Read/write support for **JSON** and **CBOR** formats. (#125)
- ✨ `io.percolator`: Added support for Comet-style N- and C-terminal modifications (#121 and #131 by @ATPs)

### Changed

- ♻️ `stats.qvalues`: Set the **regular target–decoy formula** explicitly in Pyteomics when `remove_decoy=False` and apply the **+1 correction** (probability that the first excluded decoy out-scores the threshold PSM). This produces less overly conservative q-values (e.g., on `example_files/msms.txt`). (#128)
- 🏷️ **Typing**: Adopted full **MyPy** typing across the codebase. (#125)
- 👷 **CI**: Replaced file-hash–based tests for `io.idxml` with unit tests; added formatting checks. (#125)

### Removed

- 💥 Dropped support for **Python 3.9**. (#125)

### Fixed

- 🐛 `io.mzid`: Treat **MS:1001460 “unknown modification”** as a **delta-mass–designated** modification in peptidoforms so mass calculations remain possible. Previously rendered as `[unknown modification]`. (#126 by @levitsky)
- 🐛 `io.fragpipe`: Build more comprehensive **ProForma** strings using “Assigned Modifications” from FragPipe output. (fixes #123; #124 by @levitsky)
- 🐛 `peptidoform`: Allow **residue `X` with a MassModification** to indicate a **gap of known mass** per ProForma §4.2.7; resolves failures computing theoretical mass for non-natural residues (fixes #127). (#130)

## [1.4.1] - 2025-04-15

### Fixed

- Restored compatibility with older Sage versions that have no ion mobility columns (introduced in v1.4.0) (by @rodvrees in #120)

## [1.4.0] - 2025-03-06

### Added

- ✨ `io.sage`: Add parsing of ion mobility values (PR #113)

### Fixed

- 🐛 `io.percolator`: Fix bug in `PercolatorTabWriter` where style parameter was not propagated (fixes #114, PR #117)
- 📝 Docs: Explicitly set Sphinx configuration path for Read the Docs (fixes #115, PR #118)

## [1.3.0] - 2025-01-20

### Added

- ✨ `io.idxml`: Parse ion mobility from idXML files if present.
- 🐍 Added support for Python 3.12 and 3.13

### Removed

- 🐍 Removed support for Python 3.7

### Fixed

- 🐛 Fix bug introduced in #102 where dtypes were not coerced anymore by Numpy, which lead to unexpected behavior downstream (e.g., `psm_list["is_decoy"]` would return an array of objects instead of bools)
- 🩹 Fix potential downstream issues because pepxml-read PSM had `rescoring_features=None` (partially fixes #108)

## [1.2.0] - 2024-11-19

### Added

- ✨ `io.alphadia`: Read support for AlphaDIA `precursors.tsv` (#103 by @rodvrees)
- ✨ `io.fragpipe`: Read support for FragPipe `psm.tsv` (#103 by @rodvrees)
- ✨ `io.diann`: Read support for DIA-NN TSV (#103 by @rodvrees)

### Changed

- 💥 `psm_list`: When returning a PSM property across the full PSMList (e.g. `psm_list["peptidoform"]`), `np.fromiter` is now used instead of `np.array`. This fixes an issue where if all peptidoforms have the same length, a 3D array of parsed sequences (amino acids and modifications) was be returned instead of an array of `Peptidoform` object. However, this does mean that all resulting arrays will have the `object` dtype instead of the previously coerced dtypes. This might lead to issues downstream. (#102)
- ➖ `io.idxml`: Make pyOpenMS an optional dependency, working around https://github.com/OpenMS/OpenMS/issues/7600 for now. For `idxml` support, install psm_utils with the `idxml` extra dependencies. (#107 by @paretje)

### Fixed

- 🐛 `io.pepxml`: Fix modification location and mass parsing. Position had an off-by-one error and the reported mass was the sum of the residue and modification instead of the modification alone. (fixes #100, #104)

## [1.1.1] - 2024-10-01

### Fixed

- `io`: Fix Sage filename pattern for automatic file type inference
- `io.flashlfq`: Fix writing PSMs without protein accession
- `io.flashlfq`: Fix column names `Peptide Monoisotopic Mass` and `Protein Accession`.
- `io.idxml`: Fix parsing if spectra file name not present [#92](https://github.com/compomics/psm_utils/issues/92)

## [1.1.0] - 2024-09-05

### Added

- `Peptidoform`: Add `modified_sequence` property to return the modified sequence in ProForma format, but without charge state.
- `io`: Add support for reading and writing FlashLFQ generic TSV files.

## [1.0.1] - 2024-08-28

### Fixed

- `io.percolator`: Fix and improve ScanNr inferring and writing
- `io.percolator`: Infer style from file extension if not provided (enables dynamic style determination in, for instance, `convert` function).

## [1.0.0] - 2024-08-14

### Added

- Peptidoform: Allow comparison between a peptidoform and a peptidoform string; allow direct indexing with square brackets, which indexes or slices parsed_sequence (in #89)

### Fixed

- TSV: Avoid flooding logs when reading a different file format by raising exception when three consecutive rows could not be parsed (in #88)

## [0.9.1] - 2024-07-17

### Fixed

- `io.xtandem`: Fix parsing PSMs and complete protein names in XTandem (by @julianu in #83)
- `io.tsv`: Fix warning formatting when parsing TSV (by @paretje in #85)
- `io`: Fix support for mzIdentML and pepXML files from Comet (by @paretje in #87)

## [0.9.0] - 2024-05-01

### Added

- `io`: Read and write support for writing PSMs to Apache Parquet for efficient storage of PSM lists.
- `io.sage`: Support for Sage results in Parquet format (new `SageParquetReader`, renamed `SageReader` to `SageTSVReader`).

### Changed

- Upgrade Pydantic dependency to v2. The PSM `spectrum_id` field is now always coerced to a string.
- `io.proteoscape`: Use pyarrow to iteratively read from Parquet instead of first reading an entire dataframe with Pandas.
- `io.sage`: Update compatibility to Sage v0.14
- Remove temporary patch for caching Proforma modification resolvers (now in Pyteomics v4.7.2).

## [0.8.3] - 2024-04-16

### Added

- Speed up mass calculation for large datasets by caching Proforma modification resolvers.
  Temporary patch until implemented in Pyteomics (see levitsky/pyteomics#147).

## [0.8.2] - 2024-04-05

### Added

- `io.proteoscape`: Parse filename into PSM `run` field.

## [0.8.1] - 2024-03-29

### Added

- `io.proteoscape`: Allow ProteoScapeReader instantiation from Pandas DataFrame and access PSM by index.

### Fixed

- Remove accidental print statement.
- `io.idxml`: Fixed parenthesis in type hint

### Changed

- `io.idxml`: Filter OPENMS_DATA_PATH warnings (see compomics/ms2rescore#129 and OpenMS/OpenMS#7418)
- `io.proteoscape`: Rename module from TIMScore to ProteoScape.
- `io.proteoscape`: Use correct search engine score (`x_corr_score` instead of `tims_score`)

## [0.8.0] - 2024-03-27

### Added

- `io.timscore`: Add support for TIMScore Parquet files.

### Fixed

- Fixed `_csv.Error: field larger than field limit (131072)` for very large fields when reading CSV-based PSM files.
- Pinned Pyteomics version to avoid pickling issues in multithreading (to be investigated)

## [0.7.4] - 2024-03-18

### Added

- `Peptidoform`: Support adding and applying global terminal modifications. For now using a
  workaround while waiting for official support and an implementation in Pyteomics. See
  HUPO-PSI/ProForma#6.

## [0.7.3] - 2024-03-04

### Changed

- `io.xtandem`: Parse double mass modifications as double modification instead of merging and
  summing mass shifts into a single modification.
- `io.xtandem`: Avoid float formatting issues when parsing modification mass label.
- `io.xtandem`: Parse all proteins into `protein_list` instead of only the first one.
- `io.tsv`: Log error instead of raising exception when a TSV row cannot be parsed.

## [0.7.2] - 2023-11-29

### Fixed

- `io.xtandem`: Fixed bug when extracting run name (introduced in v0.7.0)

## [0.7.1] - 2023-10-30

### Added

- Tests: Added tests for \_format_number_as_string function
- Tests: Added more test cases for `peptidoform.rename_modifications` for mass modifications
- `io.xtandem`: To parse `run` value, fall back to PSM file name if run name cannot be parsed from `label` field

### Fixed

- `peptidoform.rename_modifications`: Fixed mapping of negative mass modifications
- `io.xtandem`: Fixed regular expression to parse `run` value fom XML `label` field
- `io.idxml`: Fix handling multiple types in `rescoring_features` when writing (fixes #60)

## [0.7.0] - 2023-10-25

### Added

- `io.idxml`: Write support for idXML files, including merging an existing idXML with a `PSMList` 🎉
- `io.xtandem`: New argument `score_key` to select which score to parse as `PSM.score`.
- `io.xtandem`: Parse `run` name from X!Tandem PSM files
- Docs: Add intersphinx links to other package documentation pages.

### Changed

- `io.idxml`: Use pyOpenMS instead of Pyteomics for reading idXML (~5x faster⚡)

### Fixed

- Fix reading of pepXML files without RT
- Fixed Black formatting throughout project

## [0.6.1] - 2023-10-25

### Fixed

- `io.pepxml`: Fix reading pepXML files without retention time information.

## [0.6.0] - 2023-10-19

### Added

- `io`: Added new `io.pepxml` reader

### Fixed

- Docs: Add ionbot to README.rst, fix order in API docs

## [0.5.0] - 2023-09-20

### Added

- `Peptidoform`: Added support for `iter()` and `len()` methods
- `Peptidoform`: Added support for initialization from a `pyteomics.proforma.ProForma` object
- `PSM`: Add `precursor_mz_error` property
- `PSMList`: Added support for `append()` and `extend()` methods.
- `io`: Added new `io.ionbot` reader
- `io`: Added support for Proteome Discoverer MSF files
- `io.mzid`: Parse inverse reduced ion mobility from mzid files (e.g. from PEAKS)
- `io.mzid`: Add support for user to define custom score key
- `io.mzid`: Add `Proteome Discoverer Delta Score` to known scores (with spaces, no colons)
- `io.mzid`: Allow inconsistent presence of score in PSMs in a single mzid file

### Changed

- `PSM`: Values of the `rescoring_features` dictionary are now coerced to floats
- io: Raise `PSMUtilsIOException` when passed filetype is not known
- `io`: Make io reader `read_file` method inheritable (code cleanup)
- `io.mzid`: Throw warning when no known score can be parsed from mzid file instead of error
- `io.mzid`: Move spectrum level parsing of rt and ion mobility to function
- `io.mzid`: Give `PeptideShaker PSM score` priority over other potential search engine scores (required for correct PeptideShaker mzid parsing)
- `io.percolator`: Add option to write `PSMScore` and `ChargeN` as features to PIN file. Default is now `False`.
- Formatting: Increase max line length to 99 (code formatting)

### Fixed

- `PSMList`: Fix issue where `psm_list["protein_list"]` resulted in a Numpy error due to the inconsistent shape of the lists.
- `io.tsv`: Throw more descriptive `PSMUtilsIOException` when handeling tsv errors
- `io.msamanda`: Fix support for N/C-terminal modifications
- `io.Percolator.PercolatorTabWriter`: Allow rescoring features that are not in `feature_names` (`extrasaction` is now specified in `DictWriter`)
- Use raw strings for escape characters where needed
- Fix compatibility with sqlalchemy 2.0 (move of `declarative_base`)
- online: Remove useless == True
- docs: Set newer `build>os` configuration for readthedocs.org
- CI: Upgrade Github action versions

## [0.4.1] - 2023-07-06

### Fixed

- `PSMList`: Revert comparison operator change from v0.4.0 that results in broken `calculate_qvalues()` method (E711; Numpy array, not singleton)

## [0.4.0] - 2023-07-06

### Added

- Add + operator support for `PSMList`
- Add utility functions for m/z-mass conversion in new module `psm_utils.utils`
- `peptidoform`: Catch `ProFormaError` and reraise `PeptidoformException` with invalid peptidoform in message

### Changed

- `io.msamanda`: Changed `REQUIRED_COLUMNS` to include new features from the MS Amanda output CSV file
- `io.peptide_record`Catch the `IndexError` when a modification has a position that is out of range for the peptide, and raise an `InvalidPeprecModificationError` instead.
- Rename optional dependency `doc` to `docs`
- Implement "raise from e" when applicable throughout package

### Fixed

- Added missing `io.msamanda` API docs

## [0.3.1] - 2023-06-19

### Changed

- `io.sage`: Change `spectrum_fdr` to `spectrum_q` (crf. lazear/sage#64).

## [0.3.0] - 2023-06-08

### Added

- Add reader for [Sage](https://github.com/lazear/sage) PSM files.
- `io.mzid`: Add reading/writing of PEP and q-values

### Changed

- `psm`: The default values of `PSM.provenance_data`, `PSM.metadata` and `PSM.rescoring_features` are now `dict()` instead of `None`.
- `PSMList`: Also allow Numpy integers for indexing a single PSM
- `io.mzid.MzidReader`: Attempt to parse `retention time` or `scan start time` cvParams from both SpectrumIdentificationResult as SpectrumIdentificationItem levels. Note that according to the mzIdentML specification document (v1.1.1) neither cvParams are expected to be present at either level.
- `io.mzid.MzidReader`: Prefer `spectrum title` cvParam over `spectrumID` attribute for `PSM.spectrum_id` as these titles always match to the peak list files. In this case, `spectrumID` is saved in `metadata["mzid_spectrum_id"]`. Fall back to `spectrumID` if `spectrum title` is absent.
- `io.mzid.MzidWriter`: `PSM.retention_time` is now written as cvParam `retention time` instead of `scan start time`, and to the `SpectrumIdentificationItem` level instead of the `SpectrumIdentificationResult` level, as theoretically in psm_utils, multiple PSMs for the same spectrum can have different values for `retention_time`.
- `io.mzid.MzidWriter`: Write PSM score as cvParam `search engine specific score` instead of userParam `score`.
- `io.percolator.PercolatorTabWriter`: For PIN-style files: Use `SpecId` instead of `PSMId` and write `PSMScore` and `ChargeN` columns by default.
- Filter warnings from `psims.mzmlb` on import, as `mzmlb` is not used

### Fixed

- `psm`: Fix missing qvalue and pep in docstring
- `peptidoform`: ProForma mass modifications are now correctly parsed within the `rename_modifications` function.
- `io.maxquant.MSMSReader`: Correctly parse empty `Proteins` column to `None`
- `io.percolator.PercolatorTabReader`: Correctly parse Percolator peptidoform notation if no leading or trailing amino acids are present (e.g. `.ACDK.` instead of `K.ACDK.E`).
- `io.percolator.PercolatorTabWriter`: ScanNr is now correctly written as an integer counting from the first PSM in the file.
- `io.percolator.PercolatorTabWriter`: If no protein information is present, write the peptidoform preceded by `PEP_` to the Proteins column.
- `io.idxml`: Read metadata as strings
- `io.mzid.MzidReader`: Set `PSM.retention_time` to `None` instead of `float('nan')` if missing from the PSM file.
- `io.mzid`: Fix reading of file if charge is missing
- `io.mzid`: Fix writing if protein_list is None
- `io.mzid`: Consider all `PeptideEvidence` entries for a `SpectrumIdentificationItem` to determine `is_decoy`
- `io.mzid`: Fix handling of mzIdentML files when `is_decoy` field is not present (fixes #30)
- `io.tsv`: Raise `PSMUtilsIOException` with clear error message when TSV `protein_list` cannot be read

## [0.2.3] - 2023-03-08

### Fixed

- Fix bug in `io._base_classes` (introduced in v0.2.2)
- Fix bug in TSVReader for reading TSV files with empty protein_list

## [0.2.2] - 2023-03-08

### Fixed

- `io.peptide_record`: Fix bug where provenance item `filename` was not a string
- Various minor fixes after linting

## [0.2.1] - 2023-01-17

### Added

- `Peptidoform`: Add `is_modified` property

### Fixed

- `io.mzid`: Fix issues when parsing Comet or MSAmanda-generated mzIdentML files and certain fields are missing.

## [0.2.0] - 2022-11-12

### Added

- `PSM`: Add `ion_mobility` field
- `PSMList`: Allow slicing with bool arrays (e.g., `psm_df[psm_df["retention_time"] < 2000]`)
- `rename_modifications`: Add support for fixed modifications
- Add example files
- Online: Add support for GZipped files
- Online: Add support for logarithmic score (e.g. e-values)
- Docs: Extend contributing with example contributions
- Docs: Add notes to `PSM.get_usi()` method
- Docs: Extend quickstart on PSMList
- Docs: Add "psm_utils tags" for file formats, as used in high-level read/write/convert functions
- Docs: Peptide Record: add notes on unsupported modification types; add example for C-terminal modification
- Docs: More clearly document conversion to DataFrame
- Docs: Add bioconda install instructions
- Docs: Add citation for preprint
- Tests: Added tests for PSMList `set_ranks` and `get_rank1_psms` methods

### Changed

- `PSMList`: Refactor `set_ranks` and `get_rank1_psms` methods
- Update `.vscode/settings.json`
- Typing: Replace Union with OR operator `|`
- Online: Use percentiles instead of randomly sampling for PP plot
- Docs: Force TOC-tree max depth
- Tests: Expand unit tests in general

### Fixed

- `PSMList`: Truncate `__repr__` to first five entries only, avoiding crashing notebook output
- `Peptidoform`: Minor typing fix
- `add_fixed_modifications`: Allow input as dict as well as list of tuples
- `io`: Fix issue where the `NamedTemporaryFile` for `_supports_write_psm` was seen as invalid Percolator file
- `io.convert`: pass ` progressbar` argument to class, not `write_file`
- `io.mzid`: Add more supported MS-GF score names, make SpecEValue default
- `io.peptide_record`: `spec_id` is now a required column (`spectrum_id` is also required in PSM)
- `io.peptide_record`: Fix parsing of C-terminal modifications from proforma to peprec
- `io.percolator`: Fix Percolator peptide notation writing (fixes #18)
- `io.tsv`: Fix issue where `TSVReader` would not use string type for metadata
- `io.xtandem`: Fix issue where optional arguments were not accepted by `XTandemReader`
- `io.xtandem`: Do not split spectrum title on space
- `io.xtandem`: Fix issue where optional arguments were not accepted by `XTandemReader`
- Online: Fix pi-0 diagonal calculation
- Remove obsolete to do comments in code

## [0.1.0] - 2022-10-14

### Added

- Initial version
