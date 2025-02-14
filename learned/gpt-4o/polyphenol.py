"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol contains 2 or more benzene rings, each substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for benzene ring
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')  # pattern for benzene ring
    
    # Find all benzene rings in the molecule
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if len(benzene_matches) < 2:
        return False, "Fewer than two benzene rings"
    
    # Check for hydroxyl substitutions
    rings_with_hydroxy = 0
    for benzene_ring in benzene_matches:
        # Check if there is at least one hydroxy group on this benzene ring
        for atom_idx in benzene_ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    rings_with_hydroxy += 1
                    break  # No need to check further if a phenol -OH is found
        if rings_with_hydroxy >= 2:
            return True, "Contains 2 or more benzene rings each with an OH group"
        
    return False, f"Only {rings_with_hydroxy} benzene rings have a hydroxy group"