"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a prenol composed of more than one isoprene unit, ending with a hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for isoprene unit (C(=C)C-C)
    isoprene_pattern = Chem.MolFromSmarts("C=C-C-C")
    
    # Find all isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than 1"

    # Check for connectivity of isoprene units
    visited_atoms = set()
    chain_length = 0
    for match in isoprene_matches:
        for idx in match:
            if idx not in visited_atoms:
                chain_length += 1
                visited_atoms.add(idx)
    
    if chain_length < 2 * len(isoprene_matches):
        return False, "Isoprene units not connected as a chain"

    # Identify hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    # Check if hydroxyl is terminal
    for match in hydroxyl_matches:
        oh_atom_idx = match[0]
        oh_neighbors = mol.GetAtomWithIdx(oh_atom_idx).GetNeighbors()
        if all(neighbor.GetIdx() in visited_atoms for neighbor in oh_neighbors):
            return True, "Contains more than one isoprene unit with terminal hydroxyl group"

    return False, "Hydroxyl group is not positioned correctly relative to isoprene units"