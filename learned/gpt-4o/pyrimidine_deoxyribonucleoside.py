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

    # Define a broad pattern for deoxyribose: a furan with specific hydroxyls
    deoxyribose_pattern = Chem.MolFromSmarts("C1([C@H](O)[C@@H](CO)O)OC1")
    
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"

    # Define a broader pyrimidine SMARTS to match pyrimidine rings
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncnc1|c1ccncn1|c1cnco1|c1ccnco1")
    
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine-related base found"

    # Ensure pyrimidine is correctly attached to deoxyribose
    pyrimidine_attachment_pattern = Chem.MolFromSmarts("c1cncnc1[C@H]1OC[C@@H](O1)CO|c1ccncn1[C@H]1OC[C@@H](O1)CO")
    
    if not mol.HasSubstructMatch(pyrimidine_attachment_pattern):
        return False, "Pyrimidine base not correctly connected to deoxyribose"

    return True, "Contains pyrimidine base attached to deoxyribose sugar"

# Example test case
example_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O"  # Thymidine
result, reason = is_pyrimidine_deoxyribonucleoside(example_smiles)
print(f"Result: {result}, Reason: {reason}")