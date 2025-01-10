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

    # More general steroid pattern with flexibility in rings structure
    steroid_pattern = Chem.MolFromSmarts("[CH3]1CC2(C)C3CCC4[CH2][CH2]C(C)(C)C3CCC4C2C(C)O1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid-like backbone found"

    # Check for a hydroxyl group connected to any of the rings, typically at C3
    hydroxyl_pattern = Chem.MolFromSmarts("C(O)C")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found potentially at the 3-position"

    # Additional verification to enhance accuracy could be considered here
    # such as verifying location or orientation of the hydroxyl if necessary

    return True, "Contains steroid-like backbone with a hydroxyl group"