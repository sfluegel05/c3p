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
with fusion occurring via exactly 2 atoms.
Furthermore, by “numbering” the heterocycle beginning at its heteroatom (or, if absent, at a ring carbon bearing a carbonyl group),
the two atoms that are fused to the external benzene ring should be consecutive (positions 2 and 3)
so that the next atom (position 4) bears an external aryl (benzene) substituent.
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
       start_idx: starting atom index.
       next_idx: one neighbor of start_idx that is in ring.
       
    Returns:
       ordered list of atom indices (length = size of ring).
    """
    order = [start_idx, next_idx]
    current = next_idx
    previous = start_idx
    while len(order) < len(ring):
        # among neighbors of current that are in ring and not the previous one:
        nbrs = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(current).GetNeighbors() 
                if nbr.GetIdx() in ring and nbr.GetIdx() != previous]
        if not nbrs:
            break  # should not happen in a normal ring
        next_atom = nbrs[0]
        order.append(next_atom)
        previous, current = current, next_atom
    return order

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid.
    
    Strategy:
      1. Parse the molecule.
      2. Extract ring information.
      3. Identify candidate 6-membered rings that may serve as the heterocycle.
         We accept either:
            (a) rings that have exactly one oxygen atom (the classic 1-benzopyran),
         or (b) rings that have no oxygen but show at least 2 ring atoms that have exocyclic carbonyl groups.
      4. Identify candidate benzene rings (6-membered rings in which every atom is carbon and aromatic).
      5. For every candidate heterocycle and benzene pair that are fused (sharing exactly 2 adjacent atoms and none of these being oxygen),
         “order” the heterocycle starting from its unique heteroatom (or, if absent, from a ring carbon bearing a carbonyl).
         Then check that the two shared atoms occur as consecutive neighbors (positions 2 and 3).
      6. Then verify that the next atom (i.e. the candidate “4–position”) is carbon and carries an external substituent
         that itself is either aromatic and present within a benzene ring.
      7. Return True if these conditions are found.
    
    Args:
       smiles (str): SMILES representation of the molecule.
       
    Returns:
       (bool, str): Tuple of classification result and explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    candidate_heterocycles = []  # store tuples: (set(ring), start_atom) where start_atom is either the unique oxygen (if found)
                                # or (if no oxygen) a ring atom that has an exocyclic carbonyl bond.
    candidate_benzenes = []      # rings with 6 atoms, all aromatic carbons.
    
    for ring in all_rings:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check for benzene: every atom must be carbon and aromatic.
        if all(a.GetAtomicNum() == 6 and a.GetIsAromatic() for a in atoms):
            candidate_benzenes.append(set(ring))
        # Now check if ring is a candidate heterocycle.
        # OPTION A: The “classic” case: exactly one atom in the ring is oxygen.
        o_count = sum(1 for a in atoms if a.GetAtomicNum() == 8)
        if o_count == 1:
            # choose that oxygen as the starting atom.
            for idx in ring:
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
                    candidate_heterocycles.append((set(ring), idx))
                    break
        # OPTION B: If no oxygen in the ring, we check for at least two ring atoms bearing exocyclic carbonyl groups.
        elif o_count == 0:
            carbonyl_sites = []
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Only consider carbons:
                if atom.GetAtomicNum() != 6:
                    continue
                # Check neighbors not in ring that are oxygen in a double bond.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring:
                        continue
                    # if neighbor is oxygen and bond order is 2 (i.e. a carbonyl)
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    # bond.GetBondTypeAsDouble() returns a float value for bond order if set
                    if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                        carbonyl_sites.append(idx)
                        break  # one carbonyl is enough per ring atom
            if len(set(carbonyl_sites)) >= 2:
                # In an oxidized neoflavonoid the original heterocycle oxygen is “missing.”
                # We choose as start atom one of the ring carbons that bears a carbonyl.
                start_atom = list(set(carbonyl_sites))[0]
                candidate_heterocycles.append((set(ring), start_atom))
    
    if not candidate_heterocycles:
        return False, "No 6-membered heterocycle candidate found (neither classic pyran nor oxidized variant detected)"
    if not candidate_benzenes:
        return False, "No 6-membered aromatic benzene ring found"
        
    # For each candidate heterocycle and benzene pair, check the fusion and the substituent at “4–position.”
    for hetero_ring, start_atom in candidate_heterocycles:
        for benz_ring in candidate_benzenes:
            shared = hetero_ring.intersection(benz_ring)
            if len(shared) != 2:
                continue  # must share exactly 2 atoms
            # Also reject if any shared atom is oxygen (we do not expect the fused atoms to be heteroatoms)
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in shared):
                continue
            
            # Order the heterocycle: we use start_atom (either the oxygen or a chosen carbonyl-bearing carbon)
            # and try each neighbor in the ring.
            start_atom_obj = mol.GetAtomWithIdx(start_atom)
            ring_neighbors = [nbr.GetIdx() for nbr in start_atom_obj.GetNeighbors() if nbr.GetIdx() in hetero_ring]
            if not ring_neighbors:
                continue
            for first in ring_neighbors:
                order = order_ring(mol, hetero_ring, start_atom, first)
                if len(order) != 6:
                    continue
                # In a properly ordered 6-membered ring (cyclic), the start atom (pos1) has two neighbors:
                # They are order[1] and order[-1]. We expect the two fused atoms (shared) to be those immediately
                # adjacent to the start atom – i.e. either order[1] and order[2] or (cyclically) order[-1] and order[1].
                if set(order[1:3]) == shared:
                    candidate_pos4 = order[3]
                elif set([order[-1], order[1]]) == shared:
                    candidate_pos4 = order[2]
                else:
                    continue
                
                pos4_atom = mol.GetAtomWithIdx(candidate_pos4)
                if pos4_atom.GetAtomicNum() != 6:
                    continue  # must be a carbon
                # Check that pos4_atom has an external substituent not in the heterocycle.
                has_aryl = False
                for nbr in pos4_atom.GetNeighbors():
                    if nbr.GetIdx() in hetero_ring:
                        continue
                    # The substituent atom should be aromatic
                    if not nbr.GetIsAromatic():
                        continue
                    # And ideally belong to a benzene ring: check if in some 6-member benzene ring.
                    for ring in all_rings:
                        if len(ring) == 6 and nbr.GetIdx() in ring:
                            if all(mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in ring):
                                has_aryl = True
                                break
                    if has_aryl:
                        break
                if has_aryl:
                    return True, ("Found neoflavonoid core: a 1-benzopyran (or oxidized variant) fused to a benzene ring "
                                  "with a candidate aryl substituent at the 4–position.")
    return False, "No suitable candidate 1-benzopyran (or oxidized variant) with an external aryl substituent at the 4–position was identified"

    
# Example usage for testing:
if __name__ == "__main__":
    # You may test with one of the provided SMILES strings.
    test_smiles = "CC1=CC(=CC2=C1C(CC3(O2)CC(NC(=S)N3)(C)C)C4=CC=C(C=C4)OC)O"  # 7'-hydroxy-4'-(4-methoxyphenyl)-... thione
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)