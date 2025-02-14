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
    
    # SMARTS pattern for benzene ring with an attached hydroxy group
    hydroxy_benzene_pattern = Chem.MolFromSmarts('c1c([OH])cccc1')
    
    # Find matches for benzene rings with hydroxy groups
    matches = mol.GetSubstructMatches(hydroxy_benzene_pattern)
    
    # Count unique rings from the pattern matches
    unique_rings_with_hydroxy = set()
    for match in matches:
        ring = tuple(sorted(match))
        unique_rings_with_hydroxy.add(ring)
    
    # Ensure at least 2 unique aromatic rings have an OH substitution
    if len(unique_rings_with_hydroxy) >= 2:
        return True, "Contains 2 or more benzene rings each with an OH group"
    
    return False, f"Only {len(unique_rings_with_hydroxy)} benzene rings have a hydroxy group"