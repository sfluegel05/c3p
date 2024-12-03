import pytest

from c3p.classifier import classify
from tests.conftest import PROGRAM_DIR, METAL_ATOM, LAWRENCIUM_SMILES, HYDRIDE_SMILES


@pytest.mark.parametrize("input_smiles,expected",
                         [
                             ([LAWRENCIUM_SMILES], [(METAL_ATOM, LAWRENCIUM_SMILES)]),
                             ([HYDRIDE_SMILES], []),
                             ([LAWRENCIUM_SMILES, HYDRIDE_SMILES], [(METAL_ATOM, LAWRENCIUM_SMILES)]),
                          ])
def test_classify(input_smiles, expected):
    results = list(classify(input_smiles, PROGRAM_DIR, strict=True))
    for r in results:
        print(r)
    assert len(results) == len(input_smiles)
    pos_results = [( r.class_id, r.input_smiles) for r in results if r.is_match]
    assert len(pos_results) == len(expected)
    assert pos_results == expected
