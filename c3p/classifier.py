import logging
from pathlib import Path
from typing import Iterator, Annotated, List, Union, Optional

from c3p.datamodel import ClassificationResult, SMILES_STRING
from c3p.learn import run_code
from c3p.programs import PROGRAM_DIR

logger = logging.getLogger(__name__)


def check_class_membership(smiles_list: List[SMILES_STRING], name: str, code: str, strict=False) -> ClassificationResult:
    """
    Check if the given SMILES strings belong to the class defined in the given code.

    Args:
        smiles_list:
        name:
        code:
        strict:

    Returns:

    """
    found = False
    for line in code.split("\n"):
        # use re to check if matches ^def is_(.*)\(smiles: str):
        import re
        fn_match = re.match(r"^def (is_.*)\(smiles: str\):", line)
        if fn_match:
            found = True
            function_name = fn_match.group(1)
            for smiles in smiles_list:
                try:
                    # logger.info(f"Running {function_name} in {name} for {smiles}")
                    for _, satisfies, reason, metadata in run_code(code, function_name, [smiles], []):
                        cc = metadata.get("chemical_class", {})
                        if satisfies:
                            confidence = metadata.get("precision")
                        else:
                            # NPV = TN / (TN + FN)
                            tn = metadata.get("num_true_negatives", 0)
                            fn = metadata.get("num_false_negatives", 0)
                            confidence = tn / (tn + fn) if tn + fn > 0 else 0
                        yield ClassificationResult(
                            input_smiles=smiles,
                            class_id=cc.get("id", "-"),
                            class_name=cc.get("name", "-"),
                            is_match=satisfies,
                            reason=reason,
                            confidence=confidence
                        )
                except Exception as e:
                    if strict:
                        raise e
                    logger.error(f"Error running {function_name} in {name}: {e}")
            break
    if not found:
        if strict:
            raise ValueError(f"Could not find is_ function in {name}")


def classify(
        smiles: Union[SMILES_STRING, List[SMILES_STRING]],
        program_directory: Optional[Path] = None,
        chemical_classes: Optional[List[str]] = None,
        strict=False
    ) -> Iterator[ClassificationResult]:
    """
    Classify a SMILES string or list of SMILES strings using all programs in the given directory.

    Args:
        smiles: The SMILES string to classify
        program_directory: The directory containing the programs

    Returns:
        The classification result
    """
    # find all programs in path
    if program_directory is None:
        program_directory = PROGRAM_DIR
    programs = list(program_directory.glob("*.py"))
    logger.info(f"Found {len(programs)} programs in {program_directory}")
    if chemical_classes:
        logger.info(f"Filtering for classes: {chemical_classes}")
    smiles_list = [smiles] if isinstance(smiles, str) else smiles
    # load each program
    for program in programs:
        if program.name.startswith("__"):
            continue
        chemical_name = program.name.replace(".py", "")
        if chemical_classes is not None and chemical_name not in chemical_classes:
            logger.debug(f"Skipping {chemical_name} as not in inclusion list: {chemical_classes}")
            continue
        with open(program, "r") as f:
            code = f.read()
        yield from check_class_membership(smiles_list, chemical_name, code, strict=strict)


