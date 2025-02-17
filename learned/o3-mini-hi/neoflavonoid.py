"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: Neoflavonoid
Definition: Any 1-benzopyran with an aryl substituent at position 4.
A neoflavonoid is defined here as a bicyclic system consisting of:
  - a 6-membered heterocycle (the “pyran” ring) that contains exactly one oxygen,
  - fused to a 6-membered benzene ring (all aromatic carbons),
with fusion occurring via exactly 2 atoms.
Furthermore, by canonically “numbering” the heterocycle starting at its oxygen,
we require that the two atoms fused to the benzene correspond to positions 2 and 3,
so that the atom at position 4 (the next atom in the cycle) bears an external aryl (benzene)
substituent.
If these conditions are met the molecule is classified as a neoflavonoid.
"""

from rdkit import Chem

def order_ring(mol, ring, start_idx, next_idx):
    """
    Given a ring (as a set of atom indices) that we know is cyclic,
    return a unique ordering of the atoms starting at start_idx and then next_idx.
    We follow the ring until it cycles back.
    
    Args:
       mol: RDKit molecule.
       ring: set of atom indices (of the ring).
       start_idx: starting atom index (e.g. the unique oxygen).
       next_idx: one of the two neighbors of start_idx that is in ring.
       
    Returns:
       ordered list of atom indices (length = size of ring).
    """
    order = [start_idx, next_idx]
    current = next_idx
    previous = start_idx
    while len(order) < len(ring):
        # Among neighbors of current in the ring, choose one that is not the previous.
        nbrs = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(current).GetNeighbors() if nbr.GetIdx() in ring and nbr.GetIdx() != previous]
        if not nbrs:
            break  # Should not happen for a proper ring.
        next_atom = nbrs[0]
        order.append(next_atom)
        previous, current = current, next_atom
    return order

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid.
    
    Strategy:
      1. Parse the molecule and extract its ring information.
      2. Identify candidate 6-membered rings:
           - Candidate heterocycles: rings with exactly 6 atoms that contain exactly one oxygen.
           - Candidate benzene rings: rings with 6 atoms in which every atom is carbon and aromatic.
      3. For each candidate heterocycle and benzene pair, verify that they are fused (share exactly 2 atoms).
      4. For a given candidate heterocycle, attempt to “order” the cycle starting at its unique oxygen.
         Then require that the two fused atoms are consecutive neighbors after oxygen (positions 2 and 3).
      5. The next atom (position 4 in the ordered heterocycle) should have an external substituent
         that is (or is part of) a candidate benzene ring.
      6. If found, classify as neoflavonoid.
    
    Args:
       smiles (str): SMILES representation of the molecule.
       
    Returns:
       (bool, str): Tuple of classification result and explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    candidate_heterocycles = []  # rings meeting heterocycle criteria (6 atoms, exactly one oxygen)
    candidate_benzenes = []      # rings where every atom is carbon and aromatic.
    
    for ring in all_rings:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check for benzene: all atoms must be carbon and aromatic.
        if all(a.GetAtomicNum() == 6 and a.GetIsAromatic() for a in atoms):
            candidate_benzenes.append(set(ring))
        # Candidate heterocycle: must have exactly one oxygen (atomic num = 8) and at least one carbon.
        o_count = sum(1 for a in atoms if a.GetAtomicNum() == 8)
        c_count = sum(1 for a in atoms if a.GetAtomicNum() == 6)
        if o_count == 1 and c_count >= 1:
            candidate_heterocycles.append(set(ring))
    
    if not candidate_heterocycles:
        return False, "No 6-membered heterocycle with exactly one oxygen found (chromene core candidate missing)"
    if not candidate_benzenes:
        return False, "No 6-membered aromatic benzene ring found"
        
    # For each candidate pair (heterocycle, benzene) that are fused, check further conditions.
    for hetero in candidate_heterocycles:
        for benz in candidate_benzenes:
            shared = hetero.intersection(benz)
            if len(shared) != 2:
                continue  # Must share exactly 2 atoms.
            # Also, ensure that none of the shared atoms is oxygen.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in shared):
                continue
            
            # At this point, we assume hetero is the chromene core.
            # Find the unique oxygen in the heterocycle.
            oxygen_idx = None
            for idx in hetero:
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
                    oxygen_idx = idx
                    break
            if oxygen_idx is None:
                continue  # Should not happen.
            o_atom = mol.GetAtomWithIdx(oxygen_idx)
            # Get the neighbors of oxygen that are in the heterocycle.
            o_nbrs = [nbr.GetIdx() for nbr in o_atom.GetNeighbors() if nbr.GetIdx() in hetero]
            if len(o_nbrs) < 1:
                continue
            # For each of the two possible orders (starting at oxygen, going via one of its neighbors),
            # we check if the two fused atoms (shared with benzene) occur as the first neighbor and the next atom.
            for first in o_nbrs:
                order = order_ring(mol, hetero, oxygen_idx, first)
                if len(order) != 6:
                    continue
                # In a proper ordering of 6 atoms (cyclic), the oxygen (pos1) has two neighbors:
                # They are order[1] and order[-1]. We want the ordering where the two fused atoms (shared)
                # occur consecutively right after oxygen. Here we check: do order[1] and order[2] equal the shared set?
                if set(order[1:3]) != shared:
                    # Try the reverse order: if order[-1] and order[1] (cyclically adjacent to oxygen) are the shared set.
                    if set([order[-1], order[1]]) != shared:
                        continue
                    else:
                        # In this ordering, consider the neighbor after the shared pair.
                        # Since the ordering is cyclical, if order[-1] and order[1] are fused,
                        # we treat order[2] as candidate position 4.
                        candidate_pos4 = order[2]
                else:
                    candidate_pos4 = order[3]
                    
                # Now candidate_pos4 is our presumed “4–position” on the pyran ring.
                pos4_atom = mol.GetAtomWithIdx(candidate_pos4)
                # We require candidate pos4 to be carbon.
                if pos4_atom.GetAtomicNum() != 6:
                    continue
                # Now examine substituents of pos4_atom that are not in the heterocycle.
                for nbr in pos4_atom.GetNeighbors():
                    if nbr.GetIdx() in hetero:
                        continue
                    # The substituent must be aromatic.
                    if not nbr.GetIsAromatic():
                        continue
                    # And we check if this neighbor is in a 6-membered aromatic benzene ring.
                    for ring in all_rings:
                        if len(ring) == 6 and nbr.GetIdx() in ring:
                            if all(mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in ring):
                                return True, ("Found 1-benzopyran core with a fused benzene (fusion on two adjacent atoms) "
                                              "and with an aryl substituent at the candidate 4–position.")
    return False, "No suitable candidate 1-benzopyran core with an external aryl substituent at the 4–position was identified"


# Example usage for testing:
if __name__ == "__main__":
    # You may test with one of the provided SMILES strings.
    test_smiles = "O1[C@]2(O)[C@@]([C@H](C3=C1C(=C(O)C4=C3O[C@@H](CC4=O)C5=CC=CC=C5)C)C6=CC=CC=C6)(C(=O)C(C(=O)C2(C)C)(C)C)[H]"  # Leucadenone B
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)