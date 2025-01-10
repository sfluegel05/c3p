"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Iridoids usually consist of a cyclopentane ring fused to a six-membered oxygen heterocycle.
    Secoiridoids are formed by cleavage of a bond in the cyclopentane ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    if not atom_rings:
        return False, "No rings detected in molecule"

    # Build a list of rings as sets
    rings = [set(ring) for ring in atom_rings]

    found_iridoid = False

    # Iterate over pairs of rings
    for i in range(len(rings)):
        ring1 = rings[i]
        # Get ring size and atoms
        ring1_size = len(ring1)
        ring1_atoms = [mol.GetAtomWithIdx(idx) for idx in ring1]

        for j in range(i+1, len(rings)):
            ring2 = rings[j]
            # Get ring size and atoms
            ring2_size = len(ring2)
            ring2_atoms = [mol.GetAtomWithIdx(idx) for idx in ring2]

            # Check if rings are fused (share bonds)
            shared_bonds = set()
            for bond in mol.GetBonds():
                idx1 = bond.GetBeginAtomIdx()
                idx2 = bond.GetEndAtomIdx()
                if idx1 in ring1 and idx2 in ring1 and idx1 in ring2 and idx2 in ring2:
                    shared_bonds.add(bond)

            if not shared_bonds:
                continue  # Rings are not fused

            # Get shared atoms between rings
            shared_atoms = ring1 & ring2

            # Check if rings share two adjacent atoms (atoms connected by a bond)
            if len(shared_atoms) >= 2:
                shared_atom_indices = list(shared_atoms)
                adjacent = False
                for idx1 in shared_atom_indices:
                    for idx2 in shared_atom_indices:
                        if idx1 == idx2:
                            continue
                        bond = mol.GetBondBetweenAtoms(idx1, idx2)
                        if bond is not None:
                            adjacent = True
                            break
                    if adjacent:
                        break
                if not adjacent:
                    continue  # Shared atoms are not adjacent

                # Identify which ring is five-membered carbocycle and which is six-membered oxygen heterocycle
                if ring1_size == 5 and ring2_size == 6:
                    five_ring = ring1
                    five_ring_atoms = ring1_atoms
                    six_ring = ring2
                    six_ring_atoms = ring2_atoms
                elif ring1_size == 6 and ring2_size == 5:
                    five_ring = ring2
                    five_ring_atoms = ring2_atoms
                    six_ring = ring1
                    six_ring_atoms = ring1_atoms
                else:
                    continue  # Not a five and six-membered ring pair

                # Check that five-membered ring contains only carbons
                if not all(atom.GetAtomicNum() == 6 for atom in five_ring_atoms):
                    continue  # Five-membered ring contains non-carbon atoms

                # Check that six-membered ring contains exactly one oxygen
                num_oxygen = sum(1 for atom in six_ring_atoms if atom.GetAtomicNum() == 8)
                if num_oxygen != 1:
                    continue  # Six-membered ring does not contain exactly one oxygen

                # Ensure the oxygen atom is part of the ring (already ensured)

                # All criteria met
                found_iridoid = True
                return True, "Molecule contains cyclopentane ring fused to six-membered oxygen heterocycle"

    # For secoiridoids, check for six-membered oxygen-containing ring attached to specific groups
    # To reduce false positives, look for six-membered oxygen heterocycle with an attached ester or carboxyl group

    found_secoiridoid = False

    for ring in rings:
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            o_in_ring = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
            if len(o_in_ring) == 1:
                # Check for attached ester or carboxyl group
                ring_atom_indices = set(ring)
                oxygen_atom = o_in_ring[0]
                oxygen_neighbors = oxygen_atom.GetNeighbors()
                for neighbor in oxygen_neighbors:
                    if neighbor.GetIdx() not in ring_atom_indices:
                        # Check if neighbor is a carbon connected to a carbonyl group
                        if neighbor.GetAtomicNum() == 6:
                            for bond in neighbor.GetBonds():
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    other_atom = bond.GetOtherAtom(neighbor)
                                    if other_atom.GetAtomicNum() == 8:
                                        # Found carbonyl group adjacent to ring oxygen
                                        found_secoiridoid = True
                                        return True, "Molecule contains six-membered oxygen heterocycle attached to carbonyl group (possible secoiridoid)"

    if found_iridoid:
        return True, "Molecule contains cyclopentane ring fused to six-membered oxygen heterocycle"
    elif found_secoiridoid:
        return True, "Molecule contains six-membered oxygen heterocycle attached to carbonyl group (possible secoiridoid)"
    else:
        return False, "Molecule does not match iridoid monoterpenoid structural features"