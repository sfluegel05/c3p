"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is a naturally occurring compound derived from prostanoic acid,
    characterized by a C20 backbone with a cyclopentane ring and two side chains
    at specific positions, one of which ends with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclopentane ring not fused to other rings
    cp_ring = Chem.MolFromSmarts("C1CCCC1")  # Simple cyclopentane ring
    ring_info = mol.GetRingInfo()
    cp_rings = [ring for ring in ring_info.AtomRings() if len(ring) == 5]
    cp_ring_atoms = []
    for ring in cp_rings:
        ring_mol = mol.GetSubstructureMatches(Chem.PathToSubmol(mol, ring))
        if ring_mol:
            # Check if the ring is fused (atoms shared with other rings)
            if not any(ring_info.NumAtomRings(atom_idx) > 1 for atom_idx in ring):
                cp_ring_atoms.append(ring)
    if not cp_ring_atoms:
        return False, "No standalone cyclopentane ring found"

    # Assume the first standalone cyclopentane ring is the prostaglandin core
    cp_ring_atom_indices = cp_ring_atoms[0]

    # Check for side chains attached at positions 2 and 5 of the cyclopentane ring
    # Positions are relative, so we check that there are two side chains attached to non-adjacent ring atoms
    side_chain_atoms = []
    for idx in cp_ring_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in cp_ring_atom_indices:
                side_chain_atoms.append((idx, neighbor.GetIdx()))
    if len(side_chain_atoms) < 2:
        return False, "Less than two side chains attached to cyclopentane ring"
    # Ensure side chains are attached to non-adjacent atoms
    side_chain_positions = [atom_idx for atom_idx, _ in side_chain_atoms]
    if abs(side_chain_positions[0] - side_chain_positions[1]) == 1 or \
       abs(side_chain_positions[0] - side_chain_positions[1]) == len(cp_ring_atom_indices) - 1:
        return False, "Side chains attached to adjacent positions on cyclopentane ring"

    # Check that one side chain ends with a carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O,H1,-1]")
    has_carboxylic_acid = False
    for _, neighbor_idx in side_chain_atoms:
        side_chain = Chem.PathToSubmol(mol, Chem.FindShortestPath(mol, neighbor_idx, neighbor_idx))
        if side_chain.HasSubstructMatch(carboxylic_acid):
            has_carboxylic_acid = True
            break
    if not has_carboxylic_acid:
        return False, "No carboxylic acid group found in side chains"

    # Check for hydroxyl groups on the cyclopentane ring
    hydroxyl_on_ring = False
    for idx in cp_ring_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0:
                hydroxyl_on_ring = True
                break
    if not hydroxyl_on_ring:
        return False, "No hydroxyl groups on cyclopentane ring found"

    # Check total number of carbons (allowing some flexibility)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 17 or num_carbons > 24:
        return False, f"Number of carbons is {num_carbons}, which is not in the expected range for prostaglandins"

    # Ensure the molecule is linear (no additional rings besides the cyclopentane ring)
    if ring_info.NumRings() > 1:
        return False, "Additional rings found in the molecule"

    # If all checks pass, it is likely a prostaglandin
    return True, "Molecule matches prostaglandin structural features"

__metadata__ = {
    'chemical_class': {
        'name': 'prostaglandin',
        'definition': 'Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'error': '',
    'stdout': None
}