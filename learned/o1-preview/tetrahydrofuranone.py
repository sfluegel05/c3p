"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is any oxolane (tetrahydrofuran ring) having an oxo substituent (C=O) 
    at any position on the ring. This includes lactones where the carbonyl is part of the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()

    # Iterate over all rings in the molecule
    for idx_ring, ring in enumerate(atom_rings):
        # Check if the ring is 5-membered
        if len(ring) != 5:
            continue

        # Initialize counters
        oxygen_count = 0
        carbon_count = 0
        other_atom = False

        # Check for ring fusion
        fused = False
        for idx in ring:
            if ring_info.NumAtomRings(idx) > 1:
                fused = True
                break
        if fused:
            continue  # Skip fused rings

        # Count atoms in the ring and check for unsaturation
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        for atom in ring_atoms:
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                carbon_count += 1
            else:
                # Found an atom other than carbon or oxygen in ring
                other_atom = True
                break
            # Check if atom is aromatic
            if atom.GetIsAromatic():
                other_atom = True
                break
        if other_atom:
            continue

        # Continue only if ring has 1 oxygen, 4 carbons
        if oxygen_count != 1 or carbon_count != 4:
            continue

        # Check for unsaturation within the ring (exclude aromatic or double bonds)
        bonds_in_ring = [mol.GetBondWithIdx(bond_idx) for bond_idx in bond_rings[idx_ring]]
        ring_unsaturation = False
        for bond in bonds_in_ring:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not (
                (bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 8) or
                (bond.GetBeginAtom().GetAtomicNum() == 8 and bond.GetEndAtom().GetAtomicNum() == 6)
            ):
                # Found double bond in ring that is not C=O
                ring_unsaturation = True
                break
            # Check if bond is aromatic
            if bond.GetIsAromatic():
                ring_unsaturation = True
                break
        if ring_unsaturation:
            continue

        # Now check for presence of C=O bond in ring (lactone) or attached to ring carbons
        found_C_O_double_in_ring = False
        for bond in bonds_in_ring:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                if ((begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or
                    (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6)):
                    found_C_O_double_in_ring = True
                    break

        # Now check for C=O attached to ring carbons (outside the ring)
        found_C_O_double_attached = False
        for atom in ring_atoms:
            if atom.GetAtomicNum() != 6:
                continue  # Skip non-carbon atoms
            for neighbor in atom.GetNeighbors():
                if neighbor in ring_atoms:
                    continue  # Skip atoms in the ring
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    found_C_O_double_attached = True
                    break
            if found_C_O_double_attached:
                break

        if found_C_O_double_in_ring or found_C_O_double_attached:
            # Match found
            return True, "Contains non-fused, saturated tetrahydrofuran ring with oxo substituent"
    # No matching ring found
    return False, "No tetrahydrofuranone ring found"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:[ID]',
        'name': 'tetrahydrofuranone',
        'definition': 'Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.'
    }
}