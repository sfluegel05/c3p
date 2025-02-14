"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: cyclohexenones
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring info
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Flag to indicate cyclohexenone found
    cyclohexenone_found = False

    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if ring is six-membered
        if len(ring) != 6:
            continue

        # Check that ring is alicyclic (non-aromatic) and composed of carbon atoms
        is_valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic():
                is_valid_ring = False
                break
            if atom.GetAtomicNum() != 6:
                is_valid_ring = False
                break
            # Check that atom is only in one ring (to avoid fused rings)
            if len(atom.GetOwningMol().GetRingInfo().AtomRingsForAtom(idx)) > 1:
                is_valid_ring = False
                break
        if not is_valid_ring:
            continue

        # Collect bonds in the ring
        bond_idxs = []
        bond_atoms = []
        for i in range(len(ring)):
            atom1_idx = ring[i]
            atom2_idx = ring[(i + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
            if bond is not None:
                bond_idxs.append(bond.GetIdx())
                bond_atoms.append((atom1_idx, atom2_idx))

        # Count double bonds in the ring
        double_bond_count = 0
        for bond_idx in bond_idxs:
            bond = mol.GetBondWithIdx(bond_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1

        if double_bond_count != 1:
            continue  # Need exactly one double bond in the ring

        # Check for ketone group attached to ring carbon (exocyclic ketone)
        ketone_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check if carbon atom is double-bonded to an oxygen atom outside the ring
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetIdx() in ring:
                    continue  # Skip atoms in the ring
                if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        ketone_found = True
                        break
            if ketone_found:
                break

        if ketone_found:
            cyclohexenone_found = True
            break  # No need to check other rings

    if cyclohexenone_found:
        return True, "Contains a six-membered carbocyclic alicyclic ring with one double bond and an exocyclic ketone group"
    else:
        return False, "Does not contain a cyclohexenone ring structure"