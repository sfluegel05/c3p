"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde using its SMILES string.
    An arenecarbaldehyde is an aldehyde where the carbonyl group is attached to an aromatic moiety.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure aromaticity is set
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    
    # Check for aldehyde group pattern (C=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    # Check for aromatic moiety directly attached to the aldehyde carbon
    aromatic_ring_pattern = Chem.MolFromSmarts("[#6]a1aaaaa1")
    
    for match in aldehyde_matches:
        aldehyde_carbon_idx = match[0]
        for neighbor in mol.GetAtomWithIdx(aldehyde_carbon_idx).GetNeighbors():
            if neighbor.GetIsAromatic() and mol.HasSubstructMatch(aromatic_ring_pattern):
                return True, "Contains aldehyde group directly attached to an aromatic moiety"
    
    return False, "Aldehyde group not properly attached to aromatic moiety"