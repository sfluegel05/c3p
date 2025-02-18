"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies chemical entities of the class 1,2,4-triazines.
Definition: Any compound with a 1,2,4-triazine skeleton – a six-membered ring that is aromatic,
composed solely of carbon and nitrogen atoms (with exactly three of each) and in which some cyclic
permutation (in either forward or reverse order) of the connectivity gives the pattern:
N, N, C, N, C, C.
Improvements over the previous attempt:
  - Only aromatic rings are considered (both atoms and bonds).
  - After ordering the ring atoms by connectivity, the cyclic bond connectivity is verified.
  - The cyclic permutations (both forward and reverse) are checked for the precise pattern.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine skeleton based on its SMILES string.
    The skeleton is defined as a six-membered aromatic ring composed solely of carbon and nitrogen atoms 
    (with exactly three of each) and for which some cyclic permutation (in either forward or reverse order)
    of the ring atoms gives the pattern: N, N, C, N, C, C.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a 1,2,4-triazine core is found, False otherwise.
        str: Reason describing the classification outcome.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    # rings is a list of tuples of atom indices (order not guaranteed)
    rings = ring_info.AtomRings()
    
    # Helper function: extract an ordered ring (list of atom objects) given a set of ring indices.
    def get_ordered_ring(mol, ring_indices):
        ring_set = set(ring_indices)
        # Build neighbor mapping for atoms in the ring (neighbors that are also in the ring)
        neigh_dict = {}
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            neigh_idxs = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring_set]
            neigh_dict[idx] = neigh_idxs
        
        # Start with an arbitrary atom (first in ring_indices) and choose one neighbor to start.
        start = ring_indices[0]
        if len(neigh_dict[start]) != 2:
            return None  # not a proper ring
        for first in neigh_dict[start]:
            ordered = [start, first]
            # Iteratively pick the next neighbor (not equal to the previous in the list)
            while True:
                current = ordered[-1]
                prev = ordered[-2]
                nbrs = neigh_dict[current]
                # choose neighbor that is not the previous
                next_atom = None
                for nb in nbrs:
                    if nb != prev:
                        next_atom = nb
                        break
                if next_atom is None:
                    break
                # if we reached back at the start, ring is complete
                if next_atom == start:
                    break
                ordered.append(next_atom)
                if len(ordered) > len(ring_indices):
                    break
            if len(ordered) == len(ring_indices):
                # return atoms in the order
                return [mol.GetAtomWithIdx(idx) for idx in ordered]
        return None

    # Function to check if an ordered ring (list of atoms) matches any cyclic permutation
    # (forward or reverse) of the pattern: [N, N, C, N, C, C] (atomic numbers: 7,7,6,7,6,6).
    def matches_pattern(ordered_atoms):
        pattern = [7, 7, 6, 7, 6, 6]  # (N=7, C=6)
        n = len(pattern)
        # Test forward rotations
        for r in range(n):
            rotated = [ordered_atoms[(r + i) % n] for i in range(n)]
            if all(rotated[i].GetAtomicNum() == pattern[i] for i in range(n)):
                return True
        # Test reversed order rotations
        reversed_atoms = list(reversed(ordered_atoms))
        for r in range(n):
            rotated = [reversed_atoms[(r + i) % n] for i in range(n)]
            if all(rotated[i].GetAtomicNum() == pattern[i] for i in range(n)):
                return True
        return False

    # Iterate over each ring candidate
    for ring in rings:
        if len(ring) != 6:
            continue  # only consider six-membered rings
        
        ordered_atoms = get_ordered_ring(mol, list(ring))
        if ordered_atoms is None or len(ordered_atoms) != 6:
            continue
        
        # Check that all atoms are either carbon or nitrogen.
        if any(atom.GetAtomicNum() not in (6, 7) for atom in ordered_atoms):
            continue
        
        # Verify that exactly three nitrogens and three carbons are present.
        n_count = sum(1 for atom in ordered_atoms if atom.GetAtomicNum() == 7)
        c_count = sum(1 for atom in ordered_atoms if atom.GetAtomicNum() == 6)
        if n_count != 3 or c_count != 3:
            continue

        # Ensure that every atom in the ring is aromatic.
        if not all(atom.GetIsAromatic() for atom in ordered_atoms):
            continue

        # Ensure that all bonds between consecutive atoms in the ordered ring are aromatic.
        all_bond_aromatic = True
        num_atoms = len(ordered_atoms)
        for i in range(num_atoms):
            a1 = ordered_atoms[i]
            a2 = ordered_atoms[(i+1) % num_atoms]
            bond = mol.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx())
            if bond is None or not bond.GetIsAromatic():
                all_bond_aromatic = False
                break
        if not all_bond_aromatic:
            continue

        # Now check the cyclic permutation for the specific N, N, C, N, C, C pattern.
        if matches_pattern(ordered_atoms):
            return True, "Found 1,2,4-triazine ring pattern"
    
    return False, "No 1,2,4-triazine ring pattern found"