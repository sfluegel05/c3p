"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside includes a pyrimidine base (uracil, thymine, cytosine)
    attached to a deoxyribose sugar structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a pyrimidine base - uracil, thymine, and cytosine separately identified
    cytosine_pattern = Chem.MolFromSmarts("N1C=CC(=NC1=O)N")
    uracil_pattern = Chem.MolFromSmarts("O=C1NC=CC(=O)N1")
    thymine_pattern = Chem.MolFromSmarts("O=C1NC=CC(=O)C1")

    if not mol.HasSubstructMatch(cytosine_pattern) and not mol.HasSubstructMatch(uracil_pattern) and not mol.HasSubstructMatch(thymine_pattern):
        return False, "No pyrimidine base (uracil, thymine, or cytosine) found"

    # SMARTS pattern for deoxyribose sugar (5-membered ring with appropriate stereochemistry)
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1C[C@H](O)[C@@H](CO)O1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No correct deoxyribose (sugar) detected"

    # Check for the connection between the pyrimidine base and the deoxyribose sugar.
    # Ensure that one nitrogen atom in the pyrimidine base is connected to the sugar.
    pyrimidine_deoxyribose_linkage = Chem.MolFromSmarts("n1cn[c|n]c1[C@H]2O[C@H]([C@H](CO)O)[C@@H]2O")
    if not mol.HasSubstructMatch(pyrimidine_deoxyribose_linkage):
        return False, "Pyrimidine base is not correctly attached to deoxyribose sugar"

    return True, "Contains pyrimidine base properly attached to deoxyribose sugar"

# Test with examples
test_smiles = [
    "Nc1ccn([C@H]2CC[C@@H](CO)O2)c(=O)n1",  # zalcitabine
    "Nc1ccn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)n1",  # 2'-deoxycytidine
]

for smiles in test_smiles:
    result, reason = is_pyrimidine_deoxyribonucleoside(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")