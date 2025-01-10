"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde group at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Hexose must have exactly 6 carbons
    if c_count != 6:
        return False, f"Found {c_count} carbons, need exactly 6"

    # Hexose must have at least 5 oxygens (common in hexoses)
    if o_count < 5:
        return False, f"Found {o_count} oxygens, need at least 5"

    # Check for aldehyde group (C=O at position 1 in linear form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CX4]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if aldehyde_matches:
        return True, "Contains aldehyde group at position 1 (aldohexose)"

    # Check for ketone group (C=O at position 2 in linear form)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    for match in ketone_matches:
        # Ensure the ketone is at position 2 in the linear form
        # This is a heuristic and may not work for all cases
        if match[0] != match[1]:
            return True, "Contains ketone group at position 2 (ketohexose)"

    # If neither aldehyde nor ketone is found, it's not a hexose
    return False, "No aldehyde or ketone group found"