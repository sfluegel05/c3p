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

    # Define SMARTS pattern for isoprene unit: allowing for double bond variations
    isoprene_pattern_options = [
        Chem.MolFromSmarts("C(=C)C"),  # Common isoprene motif
        Chem.MolFromSmarts("C(=C)CC"),
        Chem.MolFromSmarts("C=C(C)C")
    ]
    
    # Check for multiple isoprene units
    match_count = 0
    for pattern in isoprene_pattern_options:
        matches = mol.GetSubstructMatches(pattern)
        match_count += len(matches)
    
    if match_count < 2:
        return False, f"Found {match_count} isoprene units, need more than 1"

    # Check for terminal hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    # Try to ensure connectivity and positioning of the hydroxyl
    visited_atoms = set()
    chain_length = 0
    for pattern in isoprene_pattern_options:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            for idx in match:
                if idx not in visited_atoms:
                    chain_length += 1
                    visited_atoms.add(idx)
    
    # Check if chain is long enough proportional to matches
    if chain_length < match_count * 3:  # Rough estimate for chain length
        return False, "Isoprene units not connected properly in a chain"

    # Verify that at least one hydroxyl is terminal
    for oh_idx in hydroxyl_matches:
        # Check if hydroxyl is on a terminal carbon
        oh_atom_idx = oh_idx[0]
        atom = mol.GetAtomWithIdx(oh_atom_idx)
        neighbors = atom.GetNeighbors()
        carbon_count = sum(1 for atom in neighbors if atom.GetAtomicNum() == 6)
        if carbon_count == 1:
            return True, "Contains more than one isoprene unit with a terminal hydroxyl group"

    return False, "Hydroxyl group is not terminal as required"