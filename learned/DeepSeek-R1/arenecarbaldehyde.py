"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde (CHEBI:22697)
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is an aldehyde where the carbonyl group is directly attached to an aromatic ring.

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
    
    # Define the aldehyde group pattern (C=O)
    aldehyde_pattern = Chem.MolFromSmarts('[CX3]=O')
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    for match in matches:
        aldehyde_carbon_idx = match[0]
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_carbon_idx)
        
        # Get neighboring atoms excluding the oxygen
        neighbors = [neighbor for neighbor in aldehyde_carbon.GetNeighbors() if neighbor.GetSymbol() != 'O']
        
        # Aldehyde should have exactly one non-oxygen neighbor (the aromatic carbon)
        if len(neighbors) != 1:
            continue
        
        neighbor_atom = neighbors[0]
        # Check if the neighbor is part of an aromatic ring
        if neighbor_atom.GetIsAromatic():
            return True, "Aldehyde group attached to aromatic ring"
    
    return False, "No aldehyde group attached to aromatic ring"