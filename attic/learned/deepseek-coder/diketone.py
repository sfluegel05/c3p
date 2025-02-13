"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: CHEBI:50699 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is a compound that contains two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ketone pattern (C=O where the carbon is bonded to two other carbons)
    ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    
    # Find all matches for the ketone pattern
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Check if there are exactly two ketone groups
    if len(ketone_matches) == 2:
        return True, "Contains exactly two ketone functional groups"
    else:
        return False, f"Found {len(ketone_matches)} ketone groups, need exactly 2"