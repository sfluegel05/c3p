"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:35762 epoxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxide(smiles: str):
    """
    Determines if a molecule contains an epoxide group based on its SMILES string.
    An epoxide is a cyclic ether where the oxygen atom forms part of a 3-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an epoxide group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for epoxide pattern: 3-membered ring with oxygen
    # [O;R1:1]1[C;R1:2][C;R1:3]1
    # O;R1 means oxygen must be in exactly one ring
    # C;R1 means carbon must be in exactly one ring
    epoxide_pattern = Chem.MolFromSmarts("[O;R1:1]1[#6;R1:2][#6;R1:3]1")
    
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found (3-membered ring with oxygen)"
    
    # Get all matches
    matches = mol.GetSubstructMatches(epoxide_pattern)
    
    # Verify ring properties for each match
    for match in matches:
        o_idx, c1_idx, c2_idx = match
        
        # Get atoms
        o_atom = mol.GetAtomWithIdx(o_idx)
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        
        # Check that oxygen is connected to exactly two atoms
        if len(o_atom.GetNeighbors()) != 2:
            continue
            
        # Check that both carbons are connected to oxygen
        c1_neighbors = set(n.GetIdx() for n in c1_atom.GetNeighbors())
        c2_neighbors = set(n.GetIdx() for n in c2_atom.GetNeighbors())
        
        if o_idx not in c1_neighbors or o_idx not in c2_neighbors:
            continue
            
        # Verify it's a 3-membered ring by checking connectivity
        if c1_idx in c2_neighbors and o_idx in c1_neighbors and o_idx in c2_neighbors:
            return True, "Contains epoxide group (3-membered ring with oxygen)"
            
    return False, "No valid epoxide group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35762',
        'name': 'epoxide',
        'definition': 'Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.',
        'parents': ['CHEBI:33641']
    }
}