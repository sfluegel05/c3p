"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside consists of a deoxyribose sugar with a pyrimidine base.

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

    # Define deoxyribose pattern: 5-membered ring with specific OH groups, accounting for typical stereochemistry
    deoxyribose_pattern = Chem.MolFromSmarts("C1C(O)C([C@@H](CO)O)O1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose backbone found"

    # Define a generic pyrimidine pattern recognizing the core pyrimidine structure
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncnc1|c1ccncn1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found"

    # Confirm pyrimidine is connected correctly to deoxyribose
    pyrimidine_attachment_pattern = Chem.MolFromSmarts("c1cncnc1[C@H]1C(O)C([C@@H](CO)O)O1")
    if not mol.HasSubstructMatch(pyrimidine_attachment_pattern):
        return False, "Pyrimidine base not correctly connected to deoxyribose"

    return True, "Contains a deoxyribose backbone with a pyrimidine base correctly attached"

# Example test case
example_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O"  # Thymidine
result, reason = is_pyrimidine_deoxyribonucleoside(example_smiles)
print(f"Result: {result}, Reason: {reason}")