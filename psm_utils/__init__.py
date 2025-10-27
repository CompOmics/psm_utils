"""Common utilities for parsing and handling PSMs, and search engine results."""

__version__ = "1.5.0"
__all__ = ["Peptidoform", "PSM", "PSMList"]

from warnings import filterwarnings

# mzmlb is not used, so hdf5plugin is not needed
filterwarnings(
    "ignore",
    message="hdf5plugin is missing",
    category=UserWarning,
    module="psims.mzmlb",
)

from psm_utils.peptidoform import Peptidoform  # noqa: E402
from psm_utils.psm import PSM  # noqa: E402
from psm_utils.psm_list import PSMList  # noqa: E402
