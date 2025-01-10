"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General steroid-like framework: allows for typical steroid four fused rings (A/B/C/D)
    steroid_framework_pattern = Chem.MolFromSmarts("C1CC2CCC3C4CCC(=O)CC4CCC3=C2C1")
    if not mol.HasSubstructMatch(steroid_framework_pattern):
        return False, "No generalized steroid-like backbone found"
    
    # 11-oxo group: look for an oxo group attached at potential position 11 in typical frameworks
    oxo_11_specific_pattern = Chem.MolFromSmarts("C1CC2CCC3C(=O)CC4CCC3=C2C1")
    if not mol.HasSubstructMatch(oxo_11_specific_pattern):
        return False, "No oxo group found at position 11"

    return True, "Contains steroid backbone with oxo group at position 11"