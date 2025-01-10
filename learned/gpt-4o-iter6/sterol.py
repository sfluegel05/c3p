"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is defined as a 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible SMARTS pattern for steroid-like backbone: 4-ring core typical in steroids
    # This pattern captures the three 6-membered rings and one 5-membered ring common in steroids
    steroid_pattern = Chem.MolFromSmiles("C1CC2CCC3C4CCCC(C4)C3CCC2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid-like backbone found"

    # Look for a hydroxyl group [-OH] which is characteristic for sterols
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    return True, "Contains steroid-like backbone with a hydroxyl group, consistent with sterol definition"