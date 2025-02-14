"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: 3'-hydroxyflavanones
"""

from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.

    A 3'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent at position 3' of the phenyl ring (B ring) attached at position 2 of the flavanone core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the chroman-4-one core (flavanone core)
    # This pattern allows for substitutions on the rings
    chroman_4_one_smarts = 'O=C1CCOc2ccccc12'
    chroman_4_one = Chem.MolFromSmarts(chroman_4_one_smarts)
    if not mol.HasSubstructMatch(chroman_4_one):
        return False, "Chroman-4-one core not found"

    # Find matches to the chroman-4-one core
    matches = mol.GetSubstructMatches(chroman_4_one)
    if not matches:
        return False, "Chroman-4-one core not found"

    # Iterate over matches to the chroman-4-one core
    for match in matches:
        # Get the atom indices of the core
        core_atoms = set(match)

        # Identify the atom at position 2 in the chroman-4-one core
        # The order of atoms in the SMARTS pattern corresponds to the match indices
        # In the SMARTS: O=C1CCOc2ccccc12
        # Atom indices in SMARTS:
        # 0: O (ketone oxygen)
        # 1: C (carbonyl carbon, position 4)
        # 2: C (position 3)
        # 3: C (position 2)
        # 4: O (ring oxygen)
        # 5-10: Aromatic carbons of the fused A ring

        position_2_idx = match[3]  # Atom at position 2
        position_2_atom = mol.GetAtomWithIdx(position_2_idx)

        # Find neighbors of position 2 atom that are not in the core (potential B ring)
        neighbors = [atom for atom in position_2_atom.GetNeighbors() if atom.GetIdx() not in core_atoms]

        # Check if any neighbor is part of an aromatic ring (B ring)
        b_ring_atom = None
        for neighbor in neighbors:
            if neighbor.GetIsAromatic():
                b_ring_atom = neighbor
                break

        if b_ring_atom is None:
            continue  # No B ring found

        # Get the B ring atoms
        b_ring_info = mol.GetRingInfo()
        b_ring_atoms = None
        for ring in b_ring_info.AtomRings():
            if b_ring_atom.GetIdx() in ring:
                b_ring_atoms = ring
                break

        if b_ring_atoms is None or len(b_ring_atoms) != 6:
            continue  # B ring is not a six-membered ring

        # Number the B ring starting from the attachment point (position 1')
        # Get the path of atoms in the B ring starting from the attachment atom
        b_ring_atom_indices = list(b_ring_atoms)
        b_ring = Chem.PathToSubmol(mol, b_ring_atom_indices)

        # Create an atom map from the B ring indices to positions 1' to 6'
        atom_idx_to_pos = {}
        try:
            # Get shortest paths from attachment atom to other atoms in the B ring
            paths = []
            for idx in b_ring_atom_indices:
                path = Chem.rdmolops.GetShortestPath(mol, b_ring_atom.GetIdx(), idx)
                if len(path) == 1:  # Same atom
                    atom_idx_to_pos[idx] = 1  # Position 1'
                else:
                    atom_idx_to_pos[idx] = len(path)  # Positions 2', 3', etc.
                    paths.append((len(path), idx))
            # Sort the positions
            positions = sorted(atom_idx_to_pos.items(), key=lambda x: x[1])
            # Get the atom index at position 3'
            pos_3_prime_idx = next(idx for idx, pos in atom_idx_to_pos.items() if pos == 3)
        except StopIteration:
            continue  # Could not determine 3' position

        pos_3_prime_atom = mol.GetAtomWithIdx(pos_3_prime_idx)

        # Check for hydroxy group at position 3'
        has_hydroxy = False
        for neighbor in pos_3_prime_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                has_hydroxy = True
                break

        if has_hydroxy:
            return True, "Contains flavanone core with hydroxy group at 3' position of B ring"

    return False, "No hydroxy group at 3' position of B ring"

__metadata__ = {
    'chemical_class': {
        'name': "3'-hydroxyflavanones",
        'definition': "Any hydroxyflavanone with a hydroxy substituent at position 3' of the phenyl ring."
    }
}