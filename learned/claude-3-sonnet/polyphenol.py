"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol contains 2 or more benzene rings each substituted by at least one hydroxy group.

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

    # Find all aromatic rings
    ri = mol.GetRingInfo()
    if not ri.NumRings():
        return False, "No rings found"

    # SMARTS pattern for benzene ring (or aromatic ring that could be part of a larger system)
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    
    # SMARTS pattern for phenol group (aromatic carbon with OH)
    phenol_pattern = Chem.MolFromSmarts("cO")
    
    # Find all benzene rings
    benzene_rings = mol.GetSubstructMatches(benzene_pattern)
    if len(benzene_rings) < 2:
        return False, f"Found only {len(benzene_rings)} benzene rings, need at least 2"
    
    # Find all phenol groups
    phenol_groups = mol.GetSubstructMatches(phenol_pattern)
    if not phenol_groups:
        return False, "No phenol groups found"
    
    # Create sets of aromatic carbons for each ring
    ring_carbons = []
    for ring in benzene_rings:
        ring_set = set(ring)
        ring_carbons.append(ring_set)
    
    # Count rings with at least one OH group
    rings_with_oh = 0
    
    # For each ring, check if it has at least one OH group
    for ring_set in ring_carbons:
        has_oh = False
        for phenol_match in phenol_groups:
            if phenol_match[0] in ring_set:  # phenol_match[0] is the aromatic carbon
                has_oh = True
                break
        if has_oh:
            rings_with_oh += 1
    
    if rings_with_oh < 2:
        return False, f"Only {rings_with_oh} rings have OH groups, need at least 2"
    
    # Additional check to ensure aromatic system
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 12:  # minimum for two benzene rings
        return False, "Not enough aromatic atoms for two benzene rings"
        
    return True, f"Found {rings_with_oh} benzene rings with OH groups"