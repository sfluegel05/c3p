"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_arenecarbaldehyde, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Kekulize the molecule to ensure proper aromatic detection
    Chem.Kekulize(mol, clearAromaticFlags=True)
    
    # Look for aldehyde pattern: [H]C=O where carbon has exactly one H
    aldehyde_pattern = Chem.MolFromSmarts("[H][CH](=O)")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
    
    # Find all aldehyde carbons
    aldehyde_carbons = mol.GetSubstructMatches(aldehyde_pattern)
    
    # Look for aromatic rings
    # Include both carbocyclic and heterocyclic aromatic systems
    aromatic_pattern = Chem.MolFromSmarts("a:a:a:a:a:a")
    five_member_aromatic = Chem.MolFromSmarts("a:a:a:a:a")
    
    if not (mol.HasSubstructMatch(aromatic_pattern) or mol.HasSubstructMatch(five_member_aromatic)):
        return False, "No aromatic ring found"
    
    # Check if any aldehyde is directly connected to an aromatic ring
    for aldehyde_match in aldehyde_carbons:
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_match[0])
        
        # Get neighboring atoms of the aldehyde carbon
        neighbors = [neighbor for neighbor in aldehyde_carbon.GetNeighbors() 
                    if neighbor.GetAtomicNum() == 6]  # Only look at carbon neighbors
        
        for neighbor in neighbors:
            # Check if the neighbor is part of an aromatic ring
            if neighbor.GetIsAromatic():
                return True, "Contains aldehyde group directly attached to aromatic ring"
    
    return False, "Aldehyde group not directly connected to aromatic ring"