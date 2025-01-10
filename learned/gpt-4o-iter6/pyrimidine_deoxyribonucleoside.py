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
    
    # Flexible pattern for deoxyribose sugar
    deoxyribose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](CO)O[C@H]1[*]")

    # Pyrimidine base pattern recognizing substitutions (using nitrogen and carbon arrangements)
    pyrimidine_base_patterns = [
        Chem.MolFromSmarts("c1ncnc2c1ncnc2"),  # Unsubstituted pyrimidine
        Chem.MolFromSmarts("c1nc[nH]c2c1ncnc2"),  # Substituted uracil/thymine/cytosine
    ]
    
    # Check for deoxyribose sugar match
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"
    
    # Check for pyrimidine base match
    if not any(mol.HasSubstructMatch(pattern) for pattern in pyrimidine_base_patterns):
        return False, "No pyrimidine base found"

    return True, "Contains both deoxyribose sugar and a pyrimidine base"