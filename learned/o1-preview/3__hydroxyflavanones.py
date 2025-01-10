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

    # Define SMARTS pattern for the flavanone core (chroman-4-one)
    # This pattern allows for substitutions at any position
    flavanone_core_smarts = 'O=C1CCC(Oc2cccc[cH]2)C1'
    flavanone_core = Chem.MolFromSmarts(flavanone_core_smarts)
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Flavanone core not found"

    # Find matches to the flavanone core
    matches = mol.GetSubstructMatches(flavanone_core)
    if not matches:
        return False, "Flavanone core not found"

    for match in matches:
        core_atoms = list(match)
        mol_atoms = mol.GetAtoms()

        # In the flavanone core SMARTS, atom indices:
        # 0: C=O carbon (position 4)
        # 1: alpha carbon (position 3)
        # 2: beta carbon (position 2)
        # 3: Oxygen atom in ring
        # 4-9: B ring carbons (phenyl ring attached at position 2)

        # Get the atom corresponding to position 2 (beta carbon)
        position_2_idx = match[2]
        position_2_atom = mol_atoms[position_2_idx]

        # Get the B ring atoms
        b_ring_atoms = set()
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if position_2_idx in ring:
                for idx in ring:
                    if idx not in core_atoms:
                        b_ring_atoms.add(idx)

        if not b_ring_atoms:
            continue  # No B ring found

        # Number the B ring starting from the attachment point (position 1')
        b_ring_indices = list(b_ring_atoms)
        attachment_atom_idx = None
        for neighbor in position_2_atom.GetNeighbors():
            if neighbor.GetIdx() in b_ring_indices:
                attachment_atom_idx = neighbor.GetIdx()
                break
        if attachment_atom_idx is None:
            continue  # Attachment atom not found

        # Get ordered B ring atoms starting from the attachment atom
        b_ring_path = Chem.rdmolops.GetShortestPath(mol, attachment_atom_idx, attachment_atom_idx)
        b_ring_ordered = b_ring_path + [attachment_atom_idx]

        # Position 1' is the attachment point, positions increase around the ring
        # Position 3' is at index 2
        if len(b_ring_ordered) < 3:
            continue  # B ring too small

        pos_3_prime_idx = b_ring_ordered[2]
        pos_3_prime_atom = mol_atoms[pos_3_prime_idx]

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