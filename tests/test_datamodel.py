import pytest

from c3p.datamodel import CodeStatistics
from tests.conftest import TEST_PROGRAM_DIR


@pytest.mark.parametrize("name", [
    "thienopyrimidine",
    "metal_atom",
    "triglyceride",
    "thioester",
])
def test_code_statistics(name):
    with open(TEST_PROGRAM_DIR / f"{name}.py", "r") as f:
        code = f.read()
        cs = CodeStatistics.from_code(code)
        print(cs.model_dump_json(indent=2))
