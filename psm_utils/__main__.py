"""
Command line interface to psm_utils.

Convert
*******

:py:mod:`psm_utils` can convert one PSM file type into another. Simply run the ``psm_utils``
command with the ``convert`` subcommand and specify the input and output filenames. By
default, the file types will be inferred from the filename and extensions. File types
can also be specified with the ``input-filetype`` and ``output-filetype`` options.

For instance:

.. code-block:: sh

    psm_utils convert --input-filetype=maxquant --output-filetype=mzid msms.txt msms_converted.mzid


Run ``psm_utils convert --help`` for an overview of all parameters.

"""

import logging
from pathlib import Path

import click
from rich.logging import RichHandler

import psm_utils.io

logger = logging.getLogger(__name__)


def main():
    """Run the main entry point for the psm_utils CLI."""
    logging.basicConfig(
        level="NOTSET",
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_suppress=[click])],
    )

    try:
        cli()
    except Exception:
        logger.exception("psm_utils encountered an error.")


@click.group()
def cli():
    """Command line interface for psm_utils."""
    pass


@cli.command()
@click.argument(
    "input-filename",
    type=click.Path(
        exists=True,
        readable=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
)
@click.argument(
    "output-filename",
    type=click.Path(
        exists=False,
        writable=True,
        dir_okay=False,
        path_type=Path,
    ),
)
@click.option("--input-filetype", default="infer", help="Input filetype")
@click.option("--output-filetype", default="infer", help="Output filetype")
def convert(
    input_filename,
    output_filename,
    input_filetype,
    output_filetype,
):
    """Convert PSM file from one format into another."""
    logger.info("Starting file conversion...")
    psm_utils.io.convert(
        input_filename,
        output_filename,
        input_filetype=input_filetype,
        output_filetype=output_filetype,
        show_progressbar=True,
    )
    logger.info("Finished file conversion! ✅")


if __name__ == "__main__":
    main()
