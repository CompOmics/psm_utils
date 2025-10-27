#########
psm_utils
#########

Common utilities for parsing and handling peptide-spectrum matches and search
engine results in Python.


.. image:: https://img.shields.io/github/v/release/compomics/psm_utils?sort=semver&style=flat-square
   :alt: GitHub release
   :target: https://github.com/compomics/psm_utils/releases

.. image:: https://img.shields.io/pypi/v/psm-utils?style=flat-square
   :alt: PyPI
   :target: https://pypi.org/project/psm-utils/

.. image:: https://img.shields.io/conda/vn/bioconda/psm-utils?style=flat-square
   :alt: Bioconda
   :target: http://bioconda.github.io/recipes/psm-utils/README.html

.. image:: https://img.shields.io/github/actions/workflow/status/compomics/psm_utils/test.yml?branch=main&label=test&style=flat-square
   :alt: GitHub Actions tests status
   :target: https://github.com/compomics/psm_utils/actions/workflows/test.yml

.. image:: https://img.shields.io/github/actions/workflow/status/compomics/psm_utils/publish.yml?event=release&style=flat-square
   :alt: GitHub Actions build status
   :target: https://github.com/compomics/psm_utils/actions/workflows/publish.yml

.. image:: https://img.shields.io/codecov/c/github/compomics/psm_utils?style=flat-square
   :alt: Codecov
   :target: https://app.codecov.io/gh/compomics/psm_utils

.. image:: https://img.shields.io/github/license/compomics/psm_utils.svg?style=flat-square
   :alt: GitHub
   :target: https://www.apache.org/licenses/LICENSE-2.0

.. image:: https://img.shields.io/badge/follow-@compomics.com-blue?style=flat-square&logo=bluesky
   :alt: Follow @compomics.com
   :target: https://bsky.app/profile/compomics.com



About
#####

Introduction
************

psm_utils is a Python package with utilities for parsing and
handling peptide-spectrum matches (PSMs) and proteomics search engine results.
It is mainly developed to be used in Python packages developed at
`CompOmics <https://www.compomics.com>`_, such as
`MS²PIP <https://github.com/compomics/ms2pip_c>`_,
`DeepLC <https://github.com/compomics/deeplc>`_, and
`MS²Rescore <https://github.com/compomics/ms2rescore>`_,
but can be useful to anyone dealing with PSMs and PSM files. Moreover, it
provides an easy-to-use CLI and
`web server <https://psm-utils.streamlitapp.com/>`_ to
convert search engine results from
one PSM file format into another.


Goals and non-goals
*******************
- To provide an easy-to-use Python API for **handling PSMs**.
- To provide a unified Python API to the plethora of **proteomics search engine
  output formats** that are in existence.
- To follow **community standards**: psm_utils pragmatically adheres to the
  standards developed by the
  `HUPO Proteomics Standards Initiative <http://psidev.info>`_, such as
  `ProForma 2.0 <https://psidev.info/proforma>`_ , the
  `Universal Spectrum Identifier <https://psidev.info/usi>`_, and
  `mzIdentML <https://psidev.info/mzidentml>`_
- To be **open and dynamic**: psm_utils is fully open source, under the
  permissive Apache 2.0 license. New reader and writer modules can easily be
  added, and we welcome everyone to contribute to the project. See
  `Contributing <https://psm-utils.readthedocs.io/en/latest/contributing>`_
  for more information.
- **NOT to reinvent the wheel**: Instead, psm_utils heavily makes
  use of packages such as `pyteomics <http://pyteomics.readthedocs.io/>`_ and
  `psims <https://github.com/mobiusklein/psims>`_ that have existing
  functionality for reading and/or writing PSM files. ``psm_utils.io``
  provides a unified, higher level Python API build on top of these packages.


Supported file formats
**********************

===================================================================================================================== ======================== =============== =============== ==========
 File format                                                                                                           psm_utils tag            Read support    Write support   Comments
===================================================================================================================== ======================== =============== =============== ==========
 `AlphaDIA precursors TSV <https://alphadia.readthedocs.io/en/latest/quickstart.html#output-files>`_                   ``alphadia``             ✅              ❌
 `CBOR <https://psm-utils.readthedocs.io/en/stable/api/psm_utils.io#module-psm_utils.io.cbor>`_                        ``cbor``                 ✅              ✅
 `DIA-NN TSV <https://github.com/vdemichev/DiaNN#output>`_                                                             ``diann``                ✅              ❌
 `FlashLFQ generic TSV <https://github.com/smith-chem-wisc/FlashLFQ/wiki/Identification-Input-Formats>`_               ``flashlfq``             ✅              ✅
 `FragPipe PSM TSV <https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html#psmtsv/>`_                       ``fragpipe``             ✅              ❌
 `ionbot CSV <https://ionbot.cloud/>`_                                                                                 ``ionbot``               ✅              ❌
 `JSON <https://psm-utils.readthedocs.io/en/stable/api/psm_utils.io#module-psm_utils.io.json>`_                        ``json``                 ✅              ✅
 `OpenMS idXML <https://www.openms.de/>`_                                                                              ``idxml``                ✅              ✅                Requires the optional ``openms`` dependency (``pip install psm-utils[idxml]``)
 `MaxQuant msms.txt <https://www.maxquant.org/>`_                                                                      ``msms``                 ✅              ❌
 `MS Amanda CSV <https://ms.imp.ac.at/?goto=msamanda>`_                                                                ``msamanda``             ✅              ❌
 `mzIdentML <https://psidev.info/mzidentml>`_                                                                          ``mzid``                 ✅              ✅
 `Parquet <https://psm-utils.readthedocs.io/en/stable/api/psm_utils.io#module-psm_utils.io.parquet>`_                  ``parquet``              ✅              ✅
 `Peptide Record <https://psm-utils.readthedocs.io/en/stable/api/psm_utils.io/#module-psm_utils.io.peptide_record>`_   ``peprec``               ✅              ✅
 `pepXML <http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML>`_                                       ``pepxml``               ✅              ❌
 `Percolator tab <https://github.com/percolator/percolator/wiki/Interface>`_                                           ``percolator``           ✅              ✅
 `Proteome Discoverer MSF <#>`_                                                                                        ``proteome_discoverer``  ✅              ❌
 `Sage Parquet <https://github.com/lazear/sage/blob/v0.14.7/DOCS.md#interpreting-sage-output>`_                        ``sage_parquet``         ✅              ❌
 `Sage TSV <https://github.com/lazear/sage/blob/v0.14.7/DOCS.md#interpreting-sage-output>`_                            ``sage_tsv``             ✅              ❌
 `ProteoScape Parquet <#>`_                                                                                            ``proteoscape``          ✅              ❌
 `TSV <https://psm-utils.readthedocs.io/en/stable/api/psm_utils.io/#module-psm_utils.io.tsv>`_                         ``tsv``                  ✅              ✅
 `X!Tandem XML <https://www.thegpm.org/tandem/>`_                                                                      ``xtandem``              ✅              ❌
===================================================================================================================== ======================== =============== =============== ==========

Legend: ✅ Supported, ❌ Unsupported



psm_utils online
################

.. image:: https://static.streamlit.io/badges/streamlit_badge_black_white.svg
   :alt: Open in streamlit
   :target: https://psm-utils.streamlitapp.com/

`psm_utils online <https://psm-utils.streamlitapp.com/>`_
is a Streamlit-based web server built on top of the psm_utils Python package. It allows
you to easily retrieve proteomics PSM statistics for any supported PSM file type, and to
convert search engine results from one PSM file format into  another. Click the badge
above to get started!



Installation
############

.. image:: https://img.shields.io/badge/install%20with-pip-brightgreen?style=flat-square
   :alt: Install with pip
   :target: https://pypi.org/project/psm-utils/

.. code-block:: sh

    pip install psm-utils


.. image:: https://img.shields.io/badge/install%20with-bioconda-blue?style=flat-square
   :alt: Install with Bioconda
   :target: http://bioconda.github.io/recipes/psm-utils/README.html

.. code-block:: sh

    conda install -c bioconda psm-utils



Full documentation
##################

The full documentation, including a quickstart guide and Python API reference
is available on `psm_utils.readthedocs.io <https://psm-utils.readthedocs.io>`_.


Citation
########

If you use psm_utils for your research, please cite the following publication:

   | **psm_utils: A high-level Python API for parsing and handling peptide-spectrum-matches and proteomics search results.**
   | Ralf Gabriels, Arthur Declercq, Robbin Bouwmeester, Sven Degroeve, Lennart Martens.
   | Journal of Proteome Research (2022). `doi:10.1021/acs.jproteome.2c00609 <https://doi.org/10.1021/acs.jproteome.2c00609>`_
