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

    # BROAD SMARTS pattern for pyrimidine base (covers uracil, thymine, cytosine features)
    pyrimidine_base_pattern = Chem.MolFromSmarts("c1cncnc1")  # Broad pattern for pyrimidine cores

    if not mol.HasSubstructMatch(pyrimidine_base_pattern):
        return False, "No broad pyrimidine core found (covers uracil, thymine, or cytosine)"

    # SMARTS pattern for a generalized deoxyribose sugar (5-membered ring with potential C or N linkage options)
    deoxyribose_pattern = Chem.MolFromSmarts("C[C@H]1O[C@@H]([C@@H](CO1)O)O")  # Generalized sugar pattern
    
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No correct deoxyribose (sugar) detected"

    # Check for the connection of base to a sugar (oxygen or nitrogen connects to sugar)
    pyrimidine_deoxyribose_linkage = Chem.MolFromSmarts("cN1CN(C=O)O1")  # Broad connection validation
    if not mol.HasSubstructMatch(pyrimidine_deoxyribose_linkage):
        return False, "Pyrimidine base is not correctly attached to deoxyribose sugar"

    return True, "Contains pyrimidine core properly attached to deoxyribose sugar"

# Test with examples
test_smiles = [
    "Nc1ccn([C@H]2CC[C@@H](CO)O2)c(=O)n1",  # zalcitabine
    "Nc1ccn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)n1",  # 2'-deoxycytidine
]

for smiles in test_smiles:
    result, reason = is_pyrimidine_deoxyribonucleoside(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")