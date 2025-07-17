from pathlib import Path
from typing import Annotated, List, Optional

import json
import typer
import yaml

import c3p.classifier as classifier
from c3p.datamodel import SMILES_STRING

app = typer.Typer(help="CHEBI classifier.")

import logging

# Set up logger
#logger = logging.getLogger(__name__)
logger = logging.getLogger()  # Get root logger


def configure_logging(verbosity: int):
    """Configure logging based on verbosity level"""
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    else:
        level = logging.WARNING

    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )

    # Configure handler
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    # Set up logger
    logger.setLevel(level)
    logger.addHandler(handler)


def verbose_option(f):
    """Decorator to add verbose option to commands"""
    return typer.Option(
        0,
        "--verbose",
        "-v",
        count=True,
        help="Verbosity level: -v for INFO, -vv for DEBUG",
    )(f)

app = typer.Typer()


@app.command()
def classify(
        smiles_list: List[SMILES_STRING],
        program_directory: Annotated[
            Optional[Path],
            typer.Option(help="path to where programs are stored")
        ] = None,
        exclude_negative: bool = typer.Option(
            False, "--exclude-negative", "-x", help="Exclude negative examples"
        ),
        confidence_cutoff: Optional[float] = typer.Option(
          0.2, "--confidence-cutoff", "-c", help="minimum value for f1_score"
        ),
        chemical_classes: Optional[List[str]] = typer.Option(
            None, "--chemical-classes", "-N", help="chemical classes to consider"
        ),
        verbose: Annotated[int, verbose_option] = 0
) -> None:
    """
    Classify SMILES strings.

    Find all classifications for a list of SMILES strings.

    Example:

        c3p "CCCCCCCCCCCCCC(O)CCCCCC"

    Only report positive matches:

        c3p -x "CCCCCCCCCCCCCC(O)CCCCCC"

    Only report positive matches with confidence above 0.75:

        c3p -x -c 0.75 "CCCCCCCCCCCCCC(O)CCCCCC"

    Pass in multiple SMILES strings:

        c3p "CCCCCCCCCCCCCC(O)CCCCCC" "CCCCCCCCCCCCCCCCCCCCO"

    TODO: fix stderr reporting.
    You may wish to redirect output

    """
    n = 0
    configure_logging(verbose)

    logger.info(f"Starting classification for SMILES: {smiles_list}")
    logger.debug(f"Searching for programs in: {program_directory}")
    chebai_style_metrics = dict()
    for result in classifier.classify(smiles_list, program_directory=program_directory, chemical_classes=chemical_classes, strict=True):
        if result.confidence is not None and confidence_cutoff is not None and result.confidence < confidence_cutoff:
            continue
        if exclude_negative and not result.is_match:
            continue
        typer.echo(yaml.dump(result.model_dump()))
        typer.echo("---")
        chebai_style_metrics[result.class_id.split(":")[-1]] = result.chebai_style_metrics
        n += 1
    if n == 0:
        typer.echo("No results found.")
    json.dump(
        chebai_style_metrics,
        open("chebai_style_metrics.json", "w"),
        indent=2,
    )
    logger.warning("PROGRAM CLASSIFICATION IS PARTIAL")
    logger.warning("RESULTS ARE EXPECTED TO BE HIGHLY INCOMPLETE")


if __name__ == "__main__":
    app()