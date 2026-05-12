"""Abstract base classes for psm_utils.io."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator
from pathlib import Path

from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class ReaderBase(ABC):
    """Abstract base class for PSM file readers."""

    filename: Path

    def __init__(
        self,
        filename: str | Path,
        *args,
        **kwargs,
    ) -> None:
        """
        Initialize PSM file reader.

        Parameters
        ----------
        filename : str or pathlib.Path
            Path to PSM file.
        *args
            Additional positional arguments for subclasses.
        **kwargs
            Additional keyword arguments for subclasses.

        """
        super().__init__()
        self.filename = Path(filename)

    def __enter__(self) -> ReaderBase:
        """Enter context manager."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Exit context manager."""
        pass

    @abstractmethod
    def __iter__(self) -> Iterator[PSM]:
        """Iterate over the PSM file and return PSMs one-by-one."""
        raise NotImplementedError()

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        return PSMList(psm_list=[psm for psm in self.__iter__()])


class WriterBase(ABC):
    """Abstract base class for PSM file writers."""

    filename: Path

    def __init__(self, filename: str | Path, *args, **kwargs) -> None:
        """
        Initialize PSM file writer.

        Parameters
        ----------
        filename : str or pathlib.Path
            Path to output PSM file.
        *args
            Additional positional arguments for subclasses.
        **kwargs
            Additional keyword arguments for subclasses.

        """
        super().__init__()
        self.filename = Path(filename)

    def __enter__(self) -> WriterBase:
        """Enter context manager."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Exit context manager."""
        pass

    @abstractmethod
    def write_psm(self, psm: PSM) -> None:
        """Write a single PSM to the PSM file."""
        raise NotImplementedError()

    @abstractmethod
    def write_file(self, psm_list: PSMList) -> None:
        """Write an entire PSMList to the PSM file."""
        raise NotImplementedError()
