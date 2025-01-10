"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfate groups
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")
    
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate groups found"

    # Look for steroid backbone pattern
    # This is a more flexible pattern to reflect the common steroid structure: 3 six-membered rings and 1 five-membered ring
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C4CCCC4CCC3C2C1")
    
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    return True, "Contains steroid backbone with sulfate group(s)"