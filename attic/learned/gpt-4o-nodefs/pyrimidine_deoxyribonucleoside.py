"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside includes a pyrimidine base (uracil, thymine, cytosine)
    attached to a deoxyribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for pyrimidine ring (aromatic 6-membered ring with 2 nitrogen atoms)
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncncc1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine ring found"

    # Check for deoxyribose sugar (detected as furanosyl ring with specific connectivities)
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1([O@H])C[C@H]([C@@H](O1)CO)O")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar detected"
    
    return True, "Contains pyrimidine base attached to deoxyribose sugar"

# Test with examples
test_smiles = [
    "Nc1ccn([C@H]2CC[C@@H](CO)O2)c(=O)n1",  # zalcitabine
    "Nc1ccn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)n1",  # 2'-deoxycytidine
]

for smiles in test_smiles:
    result, reason = is_pyrimidine_deoxyribonucleoside(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")