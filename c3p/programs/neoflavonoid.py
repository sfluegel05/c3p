"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: Neoflavonoid
Definition: Any 1-benzopyran with an aryl substituent at position 4.
A neoflavonoid is defined here as a bicyclic system consisting of:
  - a 6-membered heterocycle (the “pyran” ring) that in its “classic” form contains exactly one oxygen,
    but which in oxidized variants may have lost the ring oxygen and instead show two exocyclic carbonyl (C=O) substituents,
  - fused to a 6-membered benzene ring (all aromatic carbons),
with fusion occurring via exactly 2 atoms. Furthermore, by “numbering” the heterocycle starting from the heteroatom
(or, if none, from a ring carbon bearing a carbonyl), the two fused atoms should be consecutive – so that the next (4–)
position bears an external aryl substituent.
"""

from rdkit import Chem

def order_ring(mol, ring, start_idx, next_idx):
    """
    Given a ring (as a set of atom indices) known to be cyclic,
    return an ordered list of the atoms starting at start_idx, next_idx and following the cycle.
    """
    order = [start_idx, next_idx]
    current = next_idx
    previous = start_idx
    while len(order) < len(ring):
        nbrs = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(current).GetNeighbors() 
                if nbr.GetIdx() in ring and nbr.GetIdx() != previous]
        if not nbrs:
            break
        next_atom = nbrs[0]
        order.append(next_atom)
        previous, current = current, next_atom
    return order

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid.
    Strategy:
      1. Parse the molecule and get all 6-membered rings.
      2. Identify:
         - Candidate heterocycles: rings that either contain exactly one oxygen (classic pyran)
           or (if no oxygen is present) have at least 2 ring carbons bearing exocyclic carbonyl groups.
         - Candidate benzene rings: aromatic 6-membered rings composed entirely of carbons.
      3. For each candidate heterocycle and benzene ring pair that share exactly 2 atoms
         (and neither fused atom is oxygen), try to order (rotate) the ring so that the two fused atoms occur consecutively.
         If found, identify the next ring position (the candidate pos4).
      4. At the candidate pos4, verify that at least one neighbor not in the heterocycle is an aromatic carbon,
         and that this neighbor is part of a benzene ring that is fully external (disjoint from the heterocycle).
    Args:
       smiles (str): SMILES representation of the molecule.
    Returns:
       (bool, str): Tuple of classification result and explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # list of tuples of atom indices

    candidate_heterocycles = []  # list of tuples: (set(ring), start_atom)
    candidate_benzenes = []      # list of sets of atom indices

    # Identify candidate benzene rings and heterocycles among all 6-membered rings.
    for ring in all_rings:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Benzene: All atoms must be carbon and aromatic.
        if all(a.GetAtomicNum() == 6 and a.GetIsAromatic() for a in atoms):
            candidate_benzenes.append(set(ring))
        # Check heterocycle candidate:
        # Option A: “Classic” pyran case: exactly one oxygen in the ring.
        o_count = sum(1 for a in atoms if a.GetAtomicNum() == 8)
        if o_count == 1:
            # Use the single oxygen as the start atom.
            for idx in ring:
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
                    candidate_heterocycles.append((set(ring), idx))
                    break
        # Option B: Oxidized variant: no oxygen in the ring but at least two ring carbons that bear an exocyclic carbonyl.
        elif o_count == 0:
            carbonyl_sites = []
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring:
                        continue
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if (nbr.GetAtomicNum() == 8 and bond is not None and 
                        bond.GetBondTypeAsDouble() == 2.0):
                        carbonyl_sites.append(idx)
                        break
            if len(set(carbonyl_sites)) >= 2:
                # Pick one of the carbonyl-bearing carbons as the start atom.
                start_atom = list(set(carbonyl_sites))[0]
                candidate_heterocycles.append((set(ring), start_atom))
                
    if not candidate_heterocycles:
        return False, "No 6-membered heterocycle candidate found (neither classic pyran nor oxidized variant detected)"
    if not candidate_benzenes:
        return False, "No 6-membered aromatic benzene ring found"
        
    # For checking external aryl substituents,
    # we define a helper that verifies that a neighbor (attached to candidate pos4)
    # is part of a benzene ring that does not share any atoms with the heterocycle.
    def has_external_benzene(neighbor, hetero_ring):
        if not (neighbor.GetIsAromatic() and neighbor.GetAtomicNum() == 6):
            return False
        # Check if neighbor belongs to any 6-membered aromatic ring that is fully external to the heterocycle.
        for benz in candidate_benzenes:
            # The candidate benzene must contain this neighbor
            if neighbor.GetIdx() in benz:
                # Ensure the benzene ring is disjoint from the heterocycle.
                if hetero_ring.intersection(benz):
                    continue
                else:
                    return True
        return False

    # For each candidate heterocycle and benzene pair, check fusion and substituent.
    for hetero_ring, start_atom in candidate_heterocycles:
        for benz_ring in candidate_benzenes:
            shared_atoms = hetero_ring.intersection(benz_ring)
            if len(shared_atoms) != 2:
                continue  # need exactly 2 fused atoms
            # Reject if any fused atom is oxygen.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in shared_atoms):
                continue

            # Try ordering the heterocycle.
            start_atom_obj = mol.GetAtomWithIdx(start_atom)
            ring_neighbors = [nbr.GetIdx() for nbr in start_atom_obj.GetNeighbors() if nbr.GetIdx() in hetero_ring]
            if not ring_neighbors:
                continue
            valid_order_found = False
            candidate_pos4 = None
            for nbr in ring_neighbors:
                order = order_ring(mol, hetero_ring, start_atom, nbr)
                if len(order) != 6:
                    continue
                # Because the ring is cyclic, we try every rotation.
                n = len(order)
                for i in range(n):
                    # Let fused pair be order[i] and order[(i+1)%n].
                    pair = {order[i], order[(i+1) % n]}
                    if pair == shared_atoms:
                        # Then candidate pos4 is the next atom in order.
                        candidate_pos4 = order[(i+2) % n]
                        valid_order_found = True
                        break
                if valid_order_found:
                    break
            if not valid_order_found:
                continue

            # Now, verify the candidate pos4 has an external substituent.
            pos4_atom = mol.GetAtomWithIdx(candidate_pos4)
            # Check that at least one neighbor (via a single bond) is an aromatic carbon belonging to an external benzene.
            for nbr in pos4_atom.GetNeighbors():
                if nbr.GetIdx() in hetero_ring:
                    continue
                if has_external_benzene(nbr, hetero_ring):
                    return True, ("Found neoflavonoid core: a 1-benzopyran (or oxidized variant) fused with a benzene ring "
                                  "and with an external aryl substituent attached at the 4–position.")
    return False, "No suitable candidate 1-benzopyran (or oxidized variant) with an external aryl substituent at the 4–position was identified"


# Example usage for testing:
if __name__ == "__main__":
    # Test one of the provided SMILES strings.
    test_smiles = "CC1=CC(=CC2=C1C(CC3(O2)CC(NC(=S)N3)(C)C)C4=CC=C(C=C4)OC)O"  # 7'-hydroxy-4'-(4-methoxyphenyl)-... thione
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)