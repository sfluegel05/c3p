"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside is a deoxyribonucleoside containing a pyrimidine base.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplified deoxyribose sugar pattern
    # Recognizing a furanose ring (pentose) and typical deoxyribose features
    deoxyribose_pattern = Chem.MolFromSmarts("[C@@H]1(O)O[C@H](CO)[C@@H](O)[C@@H]1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose backbone found"

    # More flexible pyrimidine pattern accepting common variants
    # Pyrimidine base (uracil-like, thymine-like, cytosine-like)
    pyrimidine_pattern = Chem.MolFromSmarts("c1c[nH]c(=O)[nH]1 | c1cn[c,n](=[O,N])[nH]c1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found attached to the sugar"

    # Consider additional checks here if required, e.g., verifying a specific connection point

    return True, "Contains a deoxyribose backbone with a pyrimidine base"

# Test the function with example SMILES strings
example_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O"  # Thymidine
result, reason = is_pyrimidine_deoxyribonucleoside(example_smiles)
print(f"Result: {result}, Reason: {reason}")