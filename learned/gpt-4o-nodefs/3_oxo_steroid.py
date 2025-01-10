"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid has a specific steroid backbone with a ketone group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for the steroid backbone - four fused rings
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC(C4C3)C2(C1)")

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # SMARTS pattern for a 3-oxo group in the steroid backbone
    oxo_group_pattern = Chem.MolFromSmarts("CC(=O)C")
    
    # Check for 3-oxo group
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found at position 3"
    
    return True, "Contains steroid backbone with 3-oxo group"