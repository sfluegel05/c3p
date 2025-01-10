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
    
    # Define a general steroid backbone (3 six-membered and 1 five-membered condensed rings)
    steroid_pattern = Chem.MolFromSmarts("C1CC2C3CCC4CC(C)(C)C=C4C3C2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Identify an oxo group at position 11
    # 11-oxo must be associated with the specific position, typically part of ring C
    oxo_11_pattern = Chem.MolFromSmarts("C(=O)[C@H]1CC2=C([C@@H]3CCC4C(C=CC34)C12)")
    if mol.HasSubstructMatch(oxo_11_pattern):
        return True, "Structure has both the steroid backbone and an oxo group at position 11"
    
    return False, "Oxo group not accurately identified at position 11 or missing entirely"