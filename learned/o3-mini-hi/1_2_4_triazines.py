"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies chemical entities of the class 1,2,4-triazines.
Definition: Any compound with a 1,2,4-triazine skeleton,
in which the six-membered ring (whether isolated or fused) contains exactly three nitrogen atoms
and three carbon atoms, and where some cyclic permutation (or its reverse) of the connectivity
matches the pattern [N, N, C, N, C, C] (atomic numbers: [7, 7, 6, 7, 6, 6]). Additionally,
we require that every bond in the candidate ring is aromatic.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine skeleton based on its
    SMILES string. Unlike the previous attempt, fused rings are allowed as long as they
    are aromatic and the six-membered ring (as identified via SSSR) has exactly three Ns and three Cs,
    and can be cyclically rotated (or reversed) to match the pattern [7,7,6,7,6,6].
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a 1,2,4-triazine core is found, False otherwise.
        str: Reason describing the classification outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Ensure that RDKit perceives aromaticity properly.
    Chem.SanitizeMol(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each is a tuple of atom indices
    
    # Helper function: get an ordered cycle of atoms from a ring using connectivity.
    def get_ordered_ring(mol, ring_indices):
        ring_set = set(ring_indices)
        # Build a map: atom index -> list of neighbors that are also in the ring.
        neigh_dict = {}
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            neighbors_in_ring = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring_set]
            neigh_dict[idx] = neighbors_in_ring

        # Try to construct a continuous cycle;
        # start with the first atom and pick one neighbor, then extend.
        start = ring_indices[0]
        if len(neigh_dict[start]) < 2:
            return None
        for first in neigh_dict[start]:
            ordered = [start, first]
            while True:
                current = ordered[-1]
                prev = ordered[-2]
                # pick the neighbor in the ring that is not the previous atom
                next_atom = None
                for nbr in neigh_dict[current]:
                    if nbr != prev:
                        next_atom = nbr
                        break
                if next_atom is None:
                    break
                if next_atom == start:
                    # completed the cycle
                    if len(ordered) == len(ring_indices):
                        return [mol.GetAtomWithIdx(idx) for idx in ordered]
                    else:
                        break
                if next_atom in ordered:
                    break
                ordered.append(next_atom)
                if len(ordered) > len(ring_indices):
                    break
            if len(ordered) == len(ring_indices):
                # Confirm the cycle is complete by checking last and first are neighbors.
                if start in neigh_dict[ordered[-1]]:
                    return [mol.GetAtomWithIdx(idx) for idx in ordered]
        return None

    # Helper that checks whether an ordered ring (list of atoms)
    # can be rotated (or reversed and rotated) to match [7,7,6,7,6,6].
    def matches_pattern(ordered_atoms):
        target = [7, 7, 6, 7, 6, 6]
        n = len(target)
        # Get atomic numbers of the ordered atoms.
        nums = [atom.GetAtomicNum() for atom in ordered_atoms]
        # Check all cyclic rotations (forward)
        for shift in range(n):
            rotated = [nums[(shift + i) % n] for i in range(n)]
            if rotated == target:
                return True
        # Check reversed order (cyclicly rotated)
        rev = list(reversed(nums))
        for shift in range(n):
            rotated = [rev[(shift + i) % n] for i in range(n)]
            if rotated == target:
                return True
        return False

    # Iterate over each six-membered ring
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # only consider six-membered rings

        ordered_atoms = get_ordered_ring(mol, list(ring))
        if ordered_atoms is None or len(ordered_atoms) != 6:
            continue
        
        # Check composition: exactly three nitrogen and three carbon atoms.
        atomic_nums = [atom.GetAtomicNum() for atom in ordered_atoms]
        if atomic_nums.count(7) != 3 or atomic_nums.count(6) != 3:
            continue

        # Check that every connecting bond in the ring is aromatic.
        is_ring_aromatic = True
        n = len(ordered_atoms)
        for i in range(n):
            a1 = ordered_atoms[i]
            a2 = ordered_atoms[(i+1) % n]
            bond = mol.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx())
            if bond is None or not bond.GetIsAromatic():
                is_ring_aromatic = False
                break
        if not is_ring_aromatic:
            continue

        # Finally, test if the ring atoms (in some order) match our defined pattern.
        if matches_pattern(ordered_atoms):
            return True, "Found 1,2,4-triazine ring pattern"

    return False, "No 1,2,4-triazine ring pattern found"

# Example usage (for debugging/verification):
if __name__ == '__main__':
    # A small subset of examples to test the function.
    test_examples = [
        ("C1=CC=C(C=C1)C2=CN=NC(=N2)C3=CC=CC=N3", "5-phenyl-3-(2-pyridinyl)-1,2,4-triazine"),
        ("Cc1nnc(-c2ccccc2)c(=O)n1N", "metamitron"),
        ("C1=NNC(=S)NC1=O", "5-oxo-3-sulfanylidene-2H-1,2,4-triazine-6-carboxylic acid ethyl ester")
    ]
    for smi, name in test_examples:
        result, reason = is_1_2_4_triazines(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}, Reason: {reason}\n")