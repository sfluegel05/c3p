"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is defined as a molecule containing 2 or more benzene rings,
    each substituted by at least one hydroxyl group.

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

    # Define a benzene ring pattern
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if len(benzene_matches) < 2:
        return False, f"Found {len(benzene_matches)} benzene rings, need at least 2"

    # Check for hydroxyl groups on each benzene ring
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Verify that each benzene ring has at least one hydroxyl group
    for ring in benzene_matches:
        has_hydroxyl = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:  # Hydroxyl group
                    has_hydroxyl = True
                    break
            if has_hydroxyl:
                break
        if not has_hydroxyl:
            return False, f"Benzene ring {ring} has no hydroxyl group"

    return True, "Contains at least 2 benzene rings, each with at least one hydroxyl group"