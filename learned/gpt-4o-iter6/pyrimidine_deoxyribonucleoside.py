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
    
    # Flexible deoxyribose pattern
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](CO)O[C@H]1[*]")
    
    # Pyrimidine base patterns (uracil, thymine, cytosine, placeholders for derivatives)
    pyrimidine_base_patterns = [
        Chem.MolFromSmarts("C1=CN(C=O)C(=O)N=C1"),    # Generic pyrimidine pattern
        Chem.MolFromSmarts("C1=CN(C=O)C(=O)N1"),      # Uracil/derivative
        Chem.MolFromSmarts("C1=CN(C(=O)NC1)C"),       # Thymine/derivative
        Chem.MolFromSmarts("C1=CN(C(=O)NC(=N)1)C"),   # Cytidine/derivative
        Chem.MolFromSmarts("C1=CN(C(=N)N=C1)"),       # Generic cytosine modification
    ]
    
    # Check for deoxyribose sugar match
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"
    
    # Check for pyrimidine base match
    if not any(mol.HasSubstructMatch(pattern) for pattern in pyrimidine_base_patterns):
        return False, "No pyrimidine base found"

    return True, "Contains both deoxyribose sugar and pyrimidine base"