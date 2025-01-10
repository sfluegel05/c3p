"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as any aldehyde with the carbonyl group attached to an aromatic moiety.

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

    # Look for the aldehyde group (C=O attached to something, excluding full carbonyls in esters/acids)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"

    # Iterate through matched aldehyde carbons to check for attachment to aromatic ring
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    for match in matches:
        aldehyde_carbon = match[0]
        
        # Check if aldehyde carbon is attached to aromatic
        for neighbor in mol.GetAtomWithIdx(aldehyde_carbon).GetNeighbors():
            if neighbor.GetIsAromatic():
                return True, "Aldehyde group attached to an aromatic moiety"

    return False, "Aldehyde group not attached to aromatic moiety"