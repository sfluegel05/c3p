"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside is defined as a deoxynucleoside containing a pyrimidine base.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adapted deoxyribose (2'-deoxyribose) pattern - allow flexibility in stereochemistry
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](CO)[C@H](O)[C@@H]1O")
    deoxyribose_match = mol.HasSubstructMatch(deoxyribose_pattern)
    
    # More flexible pyrimidine base pattern to include uracil, thymine, cytosine bases
    pyrimidine_base_patterns = [
        Chem.MolFromSmarts("C1=NC(=O)NC(=O)C=C1"),    # Uracil
        Chem.MolFromSmarts("C1=NC(=O)N(C)C(=O)C=C1"),  # Thymine
        Chem.MolFromSmarts("C1=NC(=O)NC(N)=C1"),      # Cytosine
        Chem.MolFromSmarts("C1=NC(=O)N(C=N)C=C1"),    # Generic pyrimidine with potential variations
    ]
    
    # Check for pyrimidine base match
    pyrimidine_base_match = any(mol.HasSubstructMatch(pattern) for pattern in pyrimidine_base_patterns)
    
    if not deoxyribose_match:
        return False, "No deoxyribose sugar found"
    
    if not pyrimidine_base_match:
        return False, "No pyrimidine base found"

    return True, "Contains both deoxyribose sugar and pyrimidine base"