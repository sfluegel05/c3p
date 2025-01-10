"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is defined as having 2 or more benzene rings each with at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define patterns for benzene rings and phenol groups (benzene with hydroxy)
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')
    phenol_pattern = Chem.MolFromSmarts('c1cc(ccc1)O')
    
    # Find all benzene rings
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    
    # Must have at least 2 benzene rings
    if len(benzene_matches) < 2:
        return False, "Fewer than 2 benzene rings found"
    
    # Collect indices of benzene rings that contain a hydroxy group
    benzene_with_hydroxy_count = set()
    for match in phenol_matches:
        benzene_with_hydroxy_count.update(set(match))
    
    # Ensure each benzene ring from benzene_matches has a hydroxy
    matching_benzenes_with_hydroxy = 0
    for benzene in benzene_matches:
        if any(idx in benzene_with_hydroxy_count for idx in benzene):
            matching_benzenes_with_hydroxy += 1
            
    # Each benzene must have one hydroxy group
    if matching_benzenes_with_hydroxy < 2:
        return False, "Not all benzene rings have a hydroxy substitution"

    return True, "Contains 2 or more benzene rings each with at least one hydroxy group"