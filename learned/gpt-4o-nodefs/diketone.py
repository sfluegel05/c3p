"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone has exactly two ketone groups (C=O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ketone pattern (C=O)
    ketone_pattern = Chem.MolFromSmarts("C=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Count the number of unique ketone groups
    unique_ketone_indices = set()
    for match in ketone_matches:
        unique_ketone_indices.update(match)

    ketone_count = len(unique_ketone_indices) // 2  # divide by 2 since each match includes two atoms
    
    if ketone_count == 2:
        return True, "Contains exactly two ketone groups"
    else:
        return False, f"Found {ketone_count} ketone groups, need exactly 2"