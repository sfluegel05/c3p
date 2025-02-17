"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies chemical entities of the class 1,2,4-triazines.
Definition: Any compound with a 1,2,4-triazine skeleton is characterized by a six-membered ring,
composed solely of carbon and nitrogen atoms (exactly three of each), and in which some cyclic permutation 
(forward or reverse) of the connectivity gives the pattern: N, N, C, N, C, C.
Additionally, we require that the 1,2,4-triazine ring is not fused with another ring.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine skeleton based on its SMILES string.
    
    The skeleton is defined as a six-membered ring having exactly three nitrogens and three carbons,
    for which some cyclic permutation (forward or reverse) of the atoms (based solely on atomic number)
    matches the pattern [7, 7, 6, 7, 6, 6] (i.e. [N, N, C, N, C, C]). In addition, only non-fused (isolated)
    rings are accepted.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a 1,2,4-triazine core is found, False otherwise.
        str: Reason describing the classification outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get general ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    bond_rings = ring_info.BondRings()  # each ring is a tuple of bond indices

    # Helper: given the set/list of atom indices for a ring,
    # return an ordered list of its atoms (in connectivity order).
    def get_ordered_ring(mol, ring_indices):
        ring_set = set(ring_indices)
        # Build a dictionary: each atom idx -> list of neighbor indices that are also in the ring.
        neigh_dict = {}
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            neighs = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring_set]
            neigh_dict[idx] = neighs

        # Try to build an ordered cycle starting with the first atom
        start = ring_indices[0]
        if len(neigh_dict[start]) != 2:
            return None  # not a proper cycle
        for first in neigh_dict[start]:
            ordered = [start, first]
            while True:
                current = ordered[-1]
                prev = ordered[-2]
                # Among neighbors of current in the ring, choose one that is not prev
                next_atom = None
                for nbr in neigh_dict[current]:
                    if nbr != prev:
                        next_atom = nbr
                        break
                if next_atom is None:
                    break
                # if we returned to the start, break out (cycle complete)
                if next_atom == start:
                    break
                ordered.append(next_atom)
                if len(ordered) > len(ring_indices):
                    break
            if len(ordered) == len(ring_indices):
                # Return the list of atom objects in this order
                return [mol.GetAtomWithIdx(idx) for idx in ordered]
        return None

    # Given an ordered list of atoms from the ring, check if some cyclic permutation (forwards or reversed)
    # matches the pattern [7, 7, 6, 7, 6, 6] (N=7, C=6).
    def matches_pattern(ordered_atoms):
        pattern = [7, 7, 6, 7, 6, 6]
        n = len(pattern)
        # Check forward rotations:
        for r in range(n):
            rotated = [ordered_atoms[(r + i) % n] for i in range(n)]
            if all(rotated[i].GetAtomicNum() == pattern[i] for i in range(n)):
                return True
        # Check reversed rotations:
        reversed_atoms = list(reversed(ordered_atoms))
        for r in range(n):
            rotated = [reversed_atoms[(r + i) % n] for i in range(n)]
            if all(rotated[i].GetAtomicNum() == pattern[i] for i in range(n)):
                return True
        return False

    # Helper to check if a ring is "isolated", i.e. not fused.
    # Here we check that none of its bonds appear in any other ring.
    def is_ring_isolated(ordered_atoms):
        n = len(ordered_atoms)
        for i in range(n):
            a1 = ordered_atoms[i]
            a2 = ordered_atoms[(i+1) % n]
            bond = mol.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx())
            if bond is None:
                return False
            bond_idx = bond.GetIdx()
            # Count in how many ring bond lists this bond index appears.
            count = sum(1 for ring in bond_rings if bond_idx in ring)
            if count > 1:
                # This bond is shared with another ring (fused)
                return False
        return True

    # Iterate over all six-membered rings in the molecule.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # only consider six-membered rings
        ordered_atoms = get_ordered_ring(mol, list(ring))
        if ordered_atoms is None or len(ordered_atoms) != 6:
            continue
        
        # Check atomic composition: exactly three nitrogens and three carbons
        atomic_nums = [atom.GetAtomicNum() for atom in ordered_atoms]
        if atomic_nums.count(7) != 3 or atomic_nums.count(6) != 3:
            continue

        # Check if the ring (ignoring aromaticity flags) matches the cyclic permutation pattern.
        if not matches_pattern(ordered_atoms):
            continue

        # Only accept rings that are isolated (non-fused)
        if not is_ring_isolated(ordered_atoms):
            continue

        return True, "Found 1,2,4-triazine ring pattern"
    
    return False, "No 1,2,4-triazine ring pattern found"

# Example usage (for debugging/verification):
if __name__ == '__main__':
    # List a few SMILES strings from the provided examples:
    test_smiles = [
        "C1=CC=C(C=C1)C2=CN=NC(=N2)C3=CC=CC=N3",  # 5-phenyl-3-(2-pyridinyl)-1,2,4-triazine (TP) -> expected True
        "CC1=NNC(=O)N(C1)\\N=C\\c1cccnc1",         # pymetrozine -> expected True
        "CN(C)C1=CC=C(C=C1)NC(=O)CN2C(=O)C3=CC4=CC=CC=C4N3C=N2"  # false positive example -> expected False
    ]
    for sm in test_smiles:
        result, reason = is_1_2_4_triazines(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")