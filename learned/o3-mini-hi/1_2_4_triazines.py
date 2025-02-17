"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies chemical entities of the class 1,2,4-triazines.
Definition: Any compound that contains a 1,2,4-triazine skeleton,
which is defined as a six-membered ring made only of carbons and nitrogens
(with exactly three of each) and that in some cyclic permutation (in either forward 
or reverse order) has the pattern: N, N, C, N, C, C.
Improvements in this version over the previous attempt:
  - We no longer require the ring to be isolated (i.e. not fused).
  - We construct the ring order (connectivity) from the bond information.
  - We test cyclic rotations in both forward and reverse orders for the specific pattern.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine skeleton based on its SMILES string.
    The skeleton is defined as a six-membered ring composed solely of carbon and nitrogen atoms
    (with exactly three of each) and such that some cyclic permutation (in either forward or reverse order)
    of the ring atoms gives the pattern: N, N, C, N, C, C.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a 1,2,4-triazine core is found, False otherwise.
        str: Reason describing the classification outcome.
    """
    # Parse SMILES into molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # Returns tuples of atom indices (unordered)
    
    # Helper function: given a set of indices (ring), return the connectivity ordered list of atoms.
    def get_ordered_ring(mol, ring_indices):
        # Build a mapping: atom_idx -> list of neighbor indices (neighbors that are in the ring)
        ring_set = set(ring_indices)
        neigh_dict = {}
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            neigh_idxs = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring_set]
            neigh_dict[idx] = neigh_idxs
        
        # Start with an arbitrary atom, then choose one neighbor to begin the path.
        start = ring_indices[0]
        # There should be exactly two neighbors in a proper ring.
        if len(neigh_dict[start]) != 2:
            return None
        first = neigh_dict[start][0]
        ordered = [start, first]
        # Iteratively pick the next neighbor that is not the previous atom.
        while True:
            current = ordered[-1]
            prev = ordered[-2]
            # Among neighbors of current, choose the one that's not prev.
            nbrs = neigh_dict[current]
            if len(nbrs) < 1:
                return None
            # Ideally there is one neighbor that hasn't been visited, aside from prev.
            if nbrs[0] == prev:
                next_atom = nbrs[1] if len(nbrs) > 1 else None
            elif len(nbrs) > 1 and nbrs[1] == prev:
                next_atom = nbrs[0]
            else:
                # If both neighbors are different from prev, pick one not already in the ordered list if possible.
                next_atom = None
                for nb in nbrs:
                    if nb != prev:
                        next_atom = nb
                        break
            if next_atom is None:
                break
            # If we reached the start, then ring is complete.
            if next_atom == start:
                break
            ordered.append(next_atom)
            if len(ordered) > len(ring_indices):
                # something went wrong
                break
        if len(ordered) != len(ring_indices):
            # Try alternative ordering: if we started with a different neighbor.
            # Reset and choose the other neighbor of start.
            start = ring_indices[0]
            first = neigh_dict[start][1] if len(neigh_dict[start])>1 else None
            if first is None:
                return None
            ordered = [start, first]
            while True:
                current = ordered[-1]
                prev = ordered[-2]
                nbrs = neigh_dict[current]
                if len(nbrs) < 1:
                    return None
                if nbrs[0] == prev:
                    next_atom = nbrs[1] if len(nbrs) > 1 else None
                elif len(nbrs) > 1 and nbrs[1] == prev:
                    next_atom = nbrs[0]
                else:
                    next_atom = None
                    for nb in nbrs:
                        if nb != prev:
                            next_atom = nb
                            break
                if next_atom is None:
                    break
                if next_atom == ordered[0]:
                    break
                ordered.append(next_atom)
                if len(ordered) > len(ring_indices):
                    break
            if len(ordered) != len(ring_indices):
                return None
        # Convert atomic indices to atom objects in order.
        return [mol.GetAtomWithIdx(idx) for idx in ordered]
    
    # Function to check if a given ordered ring (list of atoms) matches any cyclic permutation (forward or reverse)
    # of the pattern: N, N, C, N, C, C
    def matches_pattern(ordered_atoms):
        pattern = [7, 7, 6, 7, 6, 6]  # atomic numbers: N=7, C=6
        n = len(pattern)
        # Check all rotations in forward order.
        for r in range(n):
            rotated = [ordered_atoms[(r + i) % n] for i in range(n)]
            if all(rotated[i].GetAtomicNum() == pattern[i] for i in range(n)):
                return True
        # Also check all rotations in reversed order.
        reversed_atoms = list(reversed(ordered_atoms))
        for r in range(n):
            rotated = [reversed_atoms[(r + i) % n] for i in range(n)]
            if all(rotated[i].GetAtomicNum() == pattern[i] for i in range(n)):
                return True
        return False
    
    # Iterate over each ring from the molecule.
    for ring in rings:
        if len(ring) != 6:
            continue  # only consider six-membered rings
        # Get the actual ordered ring using connectivity.
        ordered_atoms = get_ordered_ring(mol, list(ring))
        if ordered_atoms is None or len(ordered_atoms) != 6:
            continue
        # Check that the ring is composed only of C and N.
        if any(atom.GetAtomicNum() not in (6, 7) for atom in ordered_atoms):
            continue
        # Verify counts.
        n_count = sum(1 for atom in ordered_atoms if atom.GetAtomicNum() == 7)
        c_count = sum(1 for atom in ordered_atoms if atom.GetAtomicNum() == 6)
        if n_count != 3 or c_count != 3:
            continue
        
        # Check the cyclic permutation for the defined pattern.
        if matches_pattern(ordered_atoms):
            return True, "Found 1,2,4-triazine ring pattern"
    
    return False, "No 1,2,4-triazine ring pattern found"