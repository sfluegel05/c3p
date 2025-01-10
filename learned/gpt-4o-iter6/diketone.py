"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is defined as a compound that contains exactly two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is considered a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define ketone group pattern accounting for different environments (aromatic, aliphatic)
    ketone_pattern = Chem.MolFromSmarts("C(=O)[C!O]")  # Carbon double-bonded to oxygen and single-bonded to non-oxygen

    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"

    # Find substructure matches for ketone groups
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Check for exactly two ketone groups, ensuring non-overlapping matches
    # Using set for collecting unique carbon atoms in ketone groups
    unique_ketone_carbons = {match[0] for match in ketone_matches}

    if len(unique_ketone_carbons) == 2:
        return True, "Contains exactly 2 ketone groups, sufficient for diketone classification"
    
    if len(unique_ketone_carbons) < 2:
        return False, f"Found {len(unique_ketone_carbons)} ketone groups, not enough for diketone"
    
    return False, f"Found {len(unique_ketone_carbons)} ketone groups, exactly 2 needed for diketone"