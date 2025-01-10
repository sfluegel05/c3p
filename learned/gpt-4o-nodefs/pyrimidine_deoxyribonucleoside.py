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

    # SMARTS pattern for pyrimidine base - a 6-membered ring with 2 nitrogen atoms
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncnc1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found"

    # SMARTS pattern for deoxyribose sugar, ensuring correct stereochemistry and connectivities.
    # This pattern looks for a furanose ring (5-membered ring) with hydroxyl configuration to mimic deoxyribose.
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1([O])[C@H]([C@H]([C@@H](CO1)O)O)CO")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No correct deoxyribose sugar detected"

    # Check for the connection between the pyrimidine base and the deoxyribose sugar.
    # Ensure that the ring nitrogen is connected to the deoxyribose.
    pyrimidine_deoxyribose_linkage = Chem.MolFromSmarts("n1cncn1[C@H]2O[C@@H]([C@H](CO)O)[C@H]2O")
    if not mol.HasSubstructMatch(pyrimidine_deoxyribose_linkage):
        return False, "Pyrimidine base not correctly attached to deoxyribose sugar"

    return True, "Contains pyrimidine base properly attached to deoxyribose sugar"

# Test with examples
test_smiles = [
    "Nc1ccn([C@H]2CC[C@@H](CO)O2)c(=O)n1",  # zalcitabine
    "Nc1ccn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)n1",  # 2'-deoxycytidine
    "CC(=O)O[C@H]1C(O)O[C@H](COP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@H]1O", # False positive ADP-D-ribose example
]

for smiles in test_smiles:
    result, reason = is_pyrimidine_deoxyribonucleoside(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")