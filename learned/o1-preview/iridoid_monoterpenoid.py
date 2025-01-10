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

    # Create list of rings with atom indices
    rings = [set(ring) for ring in atom_rings]

    # Find potential iridoid core structures
    found_iridoid = False
    for i, ring1 in enumerate(rings):
        for j, ring2 in enumerate(rings):
            if i >= j:
                continue  # Avoid duplicates and self-comparison
            shared_atoms = ring1 & ring2
            if len(shared_atoms) >= 2:
                # Check if shared atoms are connected
                shared_atom_list = list(shared_atoms)
                are_adjacent = False
                for idx1 in shared_atom_list:
                    for idx2 in shared_atom_list:
                        if idx1 == idx2:
                            continue
                        bond = mol.GetBondBetweenAtoms(idx1, idx2)
                        if bond is not None:
                            are_adjacent = True
                            break
                    if are_adjacent:
                        break
                if not are_adjacent:
                    continue  # Shared atoms are not adjacent

                # Determine sizes of rings
                ring1_size = len(ring1)
                ring2_size = len(ring2)

                # Identify which ring is five-membered and which is six-membered
                if ring1_size == 5 and ring2_size == 6:
                    five_ring = ring1
                    six_ring = ring2
                elif ring1_size == 6 and ring2_size == 5:
                    five_ring = ring2
                    six_ring = ring1
                else:
                    continue  # Need one five- and one six-membered ring

                # Check atoms in five-membered ring (should be all carbons)
                five_ring_atoms = [mol.GetAtomWithIdx(idx) for idx in five_ring]
                if not all(atom.GetAtomicNum() == 6 for atom in five_ring_atoms):
                    continue  # Five-membered ring contains non-carbon atoms

                # Check atoms in six-membered ring (should contain exactly one oxygen)
                six_ring_atoms = [mol.GetAtomWithIdx(idx) for idx in six_ring]
                num_oxygen = sum(1 for atom in six_ring_atoms if atom.GetAtomicNum() == 8)
                if num_oxygen != 1:
                    continue  # Six-membered ring does not contain exactly one oxygen

                # Ensure oxygen atom is part of ring (already done)

                # All criteria met
                found_iridoid = True
                return True, "Molecule contains cyclopentane ring fused to six-membered oxygen heterocycle"

    # For secoiridoids, check for six-membered oxygen-containing ring
    # and absence of cyclopentane ring (due to cleavage)
    # First, check if there is a six-membered ring with one oxygen
    has_six_membered_oxygen_ring = False
    for ring in rings:
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            num_oxygen = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            if num_oxygen == 1:
                has_six_membered_oxygen_ring = True
                break
    # Check for cyclopentane ring
    has_cyclopentane_ring = any(len(ring) == 5 for ring in rings)

    if has_six_membered_oxygen_ring and not found_iridoid:
        return True, "Molecule contains six-membered oxygen heterocycle (possible secoiridoid)"

    return False, "Molecule does not match iridoid monoterpenoid structural features"