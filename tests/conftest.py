from pathlib import Path

import pytest

from c3p.datamodel import EvaluationResult

THIS_DIR = Path(__file__).parent
OUTPUT_DIR = THIS_DIR / "output"
INPUT_DIR = THIS_DIR / "input"
PROGRAM_DIR = INPUT_DIR / "programs"

@pytest.fixture
def eval_result_metal_atom() -> EvaluationResult:
    with open(INPUT_DIR / "metal_atom.json", "r") as f:
        return EvaluationResult.model_validate_json(f.read())


METAL_ATOM = "CHEBI:33521"
LAWRENCIUM_SMILES = "[Lr]"
HYDRIDE_SMILES = "[H-]"
