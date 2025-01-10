"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: 3'-hydroxyflavanones
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define a more general SMARTS pattern for the flavanone core
    # Flavanone core: two fused rings (benzopyranone) with ketone at position 4
    flavanone_core_smarts = 'O=C1CCC2=CC=CC=C2O1'
    flavanone_core = Chem.MolFromSmarts(flavanone_core_smarts)
    if flavanone_core is None:
        return False, "Invalid flavanone core SMARTS"

    # Check for flavanone core
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Flavanone core not found"

    # Find the matching flavanone core in the molecule
    core_matches = mol.GetSubstructMatches(flavanone_core)
    if not core_matches:
        return False, "Flavanone core not found"

    # Iterate over core matches to find B ring and check for hydroxy at 3'
    for core_match in core_matches:
        core_atom_indices = set(core_match)

        # Find atom at position 2 (attachment point for B ring)
        # In the SMARTS pattern, atom index 3 corresponds to position 2
        position_2_idx = core_match[3]
        position_2_atom = mol.GetAtomWithIdx(position_2_idx)

        # Find neighboring atoms not in the core (B ring)
        b_ring_atom = None
        for neighbor in position_2_atom.GetNeighbors():
            if neighbor.GetIdx() not in core_atom_indices:
                b_ring_atom = neighbor
                break
        if b_ring_atom is None:
            continue  # B ring not found

        # Check if the B ring is a benzene ring
        b_ring_info = mol.GetRingInfo()
        b_ring_atoms = None
        for ring in b_ring_info.AtomRings():
            if b_ring_atom.GetIdx() in ring:
                ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                if all(atom.GetIsAromatic() for atom in ring_atoms) and len(ring_atoms) == 6:
                    b_ring_atoms = ring
                    break
        if b_ring_atoms is None:
            continue  # B ring not found

        # Number the B ring atoms starting from the attachment point
        b_ring_atom_indices = list(b_ring_atoms)
        attachment_idx_in_ring = b_ring_atom_indices.index(b_ring_atom.GetIdx())
        ordered_b_ring_atoms = b_ring_atom_indices[attachment_idx_in_ring:] + b_ring_atom_indices[:attachment_idx_in_ring]

        # Position 1' is the attachment point, positions increase around the ring
        # Position 3' is at index 2
        pos_3_prime_idx = ordered_b_ring_atoms[2]
        pos_3_prime_atom = mol.GetAtomWithIdx(pos_3_prime_idx)

        # Check for hydroxy group at position 3'
        hydroxy_group_found = False
        for neighbor in pos_3_prime_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                hydroxy_group_found = True
                break

        if hydroxy_group_found:
            return True, "Contains flavanone core with hydroxy group at 3' position of B ring"

    return False, "No hydroxy group at 3' position of B ring"

__metadata__ = {
    'chemical_class': {
        'name': "3'-hydroxyflavanones",
        'definition': "Any hydroxyflavanone with a hydroxy substituent at position 3' of the phenyl ring."
    }
}