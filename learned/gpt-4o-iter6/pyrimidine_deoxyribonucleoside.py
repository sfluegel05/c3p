"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside is defined as a deoxyribonucleoside containing a pyrimidine base.

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
    
    # More specific pattern for deoxyribose sugar, allowing stereochemistry variations
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)C[C@@H](CO)O1")

    # Revised patterns for pyrimidine bases, allowing variations and common modifications
    pyrimidine_base_patterns = [
        Chem.MolFromSmarts("c1ncnc[nH]1"),  # Generic pyrimidine ring including cytosine, thymine, uracil core
        Chem.MolFromSmarts("c1[nH]cnc1"),   # Allow hydroxyl/methylamide modified pyrimidines
        Chem.MolFromSmarts("n1c(=O)[nH]c(=O)c1") # Thymine and uracil derivatives
    ]
    
    # Check for deoxyribose sugar match
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"
    
    # Check for any pyrimidine base match
    if not any(mol.HasSubstructMatch(pattern) for pattern in pyrimidine_base_patterns):
        return False, "No pyrimidine base found"

    return True, "Contains both deoxyribose sugar and a pyrimidine base"