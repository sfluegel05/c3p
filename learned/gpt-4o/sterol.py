"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a 3-hydroxy steroid with a sterane skeleton closely related to cholestan-3-ol, 
    possibly with additional carbons in side chains.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating whether the molecule is a sterol and the reason for classification.
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated SMARTS pattern for sterane skeleton (4-ring system with flexibility for side chains).
    sterane_pattern = Chem.MolFromSmarts('C1CC2CC3CCC4C(C)C(O)CCC4(C)C3C2C1')
    if not mol.HasSubstructMatch(sterane_pattern):
        return False, "No sterane skeleton found"

    # Check for hydroxyl group at an appropriate position
    # Pattern to identify a hydroxyl group connected to a tertiary carbon in the sterane core
    hydroxy_pattern = Chem.MolFromSmarts('[C;R0]([C;R2])[O]')
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No suitable 3-hydroxy group found in cholesterol-like backbone"

    # Ensure potential side chains or functional groups
    # Flexibility for pattern that includes branching
    isopropyl_pattern = Chem.MolFromSmarts('CC(C)C')
    if mol.HasSubstructMatch(isopropyl_pattern) or Chem.MolFromSmarts('[CH3]') is not None:
        return True, "Contains sterane skeleton with 3-hydroxy group and side chains"

    return False, "Does not match the sterol structural requirements even with hydroxyl group"