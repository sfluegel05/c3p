import pytest
from typer.testing import CliRunner

from c3p.cli import app
from tests.conftest import LAWRENCIUM_SMILES, METAL_ATOM, HYDRIDE_SMILES, TEST_PROGRAM_DIR

runner = CliRunner()

@pytest.mark.parametrize("input_smiles,expected",
                        [
                            ([LAWRENCIUM_SMILES], [(METAL_ATOM, LAWRENCIUM_SMILES)]),
                            ([HYDRIDE_SMILES], []),
                            ([LAWRENCIUM_SMILES, HYDRIDE_SMILES], [(METAL_ATOM, LAWRENCIUM_SMILES)]),
                        ])
def test_classify_command(input_smiles, expected):
    result = runner.invoke(app, ["classify", *input_smiles, "--program-directory", TEST_PROGRAM_DIR])
    print(result.stdout)
    assert result.exit_code == 0
    for e in expected:
        assert e[0] in result.stdout
        assert e[1] in result.stdout