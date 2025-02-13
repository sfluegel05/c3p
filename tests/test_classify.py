import pytest

from c3p.classifier import classify
from tests.conftest import TEST_PROGRAM_DIR, METAL_ATOM, LAWRENCIUM_SMILES, HYDRIDE_SMILES

# TODO: make this less rigid; needs changes if any test programs are added
NUM_PROGRAMS = 3

@pytest.mark.parametrize("input_smiles,expected",
                         [
                             ([LAWRENCIUM_SMILES], [(METAL_ATOM, LAWRENCIUM_SMILES)]),
                             ([HYDRIDE_SMILES], []),
                             ([LAWRENCIUM_SMILES, HYDRIDE_SMILES], [(METAL_ATOM, LAWRENCIUM_SMILES)]),
                          ])
def test_classify(input_smiles, expected):
    results = list(classify(input_smiles, TEST_PROGRAM_DIR, strict=True))
    for r in results:
        print(r)
    assert len(results) == NUM_PROGRAMS * len(input_smiles)
    pos_results = [( r.class_id, r.input_smiles) for r in results if r.is_match]
    assert len(pos_results) == len(expected)
    assert pos_results == expected


