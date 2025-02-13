"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:33648 epoxide
An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxide(smiles: str):
    """
    Determines if a molecule contains an epoxide group based on its SMILES string.

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
    
    # Look for epoxide pattern ([O;r3])
    epoxide_pattern = Chem.MolFromSmarts("[O;r3]")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    
    # Check if the matches are indeed cyclic ethers
    for match in epoxide_matches:
        atom = mol.GetAtomWithIdx(match)
        if atom.GetIsAromatic() or atom.GetExplicitDegree() != 2:
            continue  # Skip non-epoxide matches
        
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            continue  # Oxygen must have exactly 2 neighbors
        
        if neighbors[0].GetAtomicNum() != 6 or neighbors[1].GetAtomicNum() != 6:
            continue  # Neighbors must be carbon atoms
        
        # Check if the carbons form a 3-membered ring with the oxygen
        ring = mol.GetAtomRingInfo().IsCyclic()
        if ring and len(ring) == 3:
            return True, "Contains a 3-membered cyclic ether (epoxide) group"
    
    return False, "No epoxide group found"