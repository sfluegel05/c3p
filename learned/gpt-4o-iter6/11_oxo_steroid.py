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
    
    # Generalized steroid framework: three cyclohexane rings and one cyclopentane ring 
    # Combined as: [6-6-6-5] with possible connections
    steroid_pattern = Chem.MolFromSmarts("C1CCC2CC3CCC4C(C3)CC2C1C4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No generalized steroid-like backbone found"
    
    # Identify an oxo group on the C ring, which is usually the third ring in a steroid structure
    oxo_11_pattern = Chem.MolFromSmarts("[C,c;R]{3}[C;R]=O") 
    if mol.HasSubstructMatch(oxo_11_pattern):
        return True, "Structure has a generalized steroid backbone with an oxo group at position 11"
    
    return False, "Oxo group not accurately identified at position 11 or missing entirely"