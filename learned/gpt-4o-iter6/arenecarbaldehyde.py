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
    
    # Define an aldehyde group pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Define a more comprehensive aromatic pattern that captures various aromatics, 
    # including benzene and heteroaromatics
    aromatic_pattern = Chem.MolFromSmarts("a")
    
    for match in aldehyde_matches:
        aldehyde_atom = match[0] # Carbon of the aldehyde
        
        # Check if connected atom is aromatic
        for neighbor in mol.GetAtomWithIdx(aldehyde_atom).GetNeighbors():
            if neighbor.HasSubstructMatch(aromatic_pattern):
                return True, "Aldehyde group attached to an aromatic moiety"
    
    # If no attachments to aromatic moiety are found:
    return False, "Aldehyde group not attached to aromatic moiety"