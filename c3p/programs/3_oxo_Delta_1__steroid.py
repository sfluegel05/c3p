"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-Delta(1) steroid.
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
Improved approach:
  1. Identify a fused steroid nucleus by analyzing the rings in the molecule.
     We require a connected set of >=4 rings (ring fusion defined by sharing at least 2 atoms)
     and that within that component at least one ring is 5-membered (as typical).
  2. Look for a ring carbonyl group ([R]C(=O)[R]) that belongs to one of the candidate rings.
  3. Instead of forcing the double bond to be directly on the carbonyl carbon, 
     we check that at least one neighbor (within the same ring) of the carbonyl carbon 
     participates in a double bond (with another atom in that ring) – an approximation of the Delta(1) feature.
     
If any of these conditions are not met we return False with the appropriate reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid.
    
    This function works in two steps:
      (i) It first identifies a fused ring system (the steroid nucleus) by collecting all ring systems and 
          then determining their connectivity. A typical steroid nucleus is a fused system with at least 
          four rings, and usually comprises one five-membered and three six-membered rings.
      (ii) It then finds a ring carbonyl group (3-oxo) within that fused system and checks if, in the same ring,
           there is at least one double bond connected to a neighbor of the carbonyl carbon, as a proxy for the 
           Delta(1) (i.e., 1,2-double bond) feature.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-oxo-Delta(1) steroid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all ring atom index lists.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # Build a graph connecting rings that share at least 2 atoms (i.e. fused rings).
    n_rings = len(rings)
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(rings[i])
        for j in range(i+1, n_rings):
            set_j = set(rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
                
    # Find connected components in the ring graph.
    visited = set()
    components = []
    for i in range(n_rings):
        if i in visited:
            continue
        comp = set()
        stack = [i]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            stack.extend(ring_graph[node] - visited)
        components.append(comp)
    
    # Look for a candidate fused ring component that resembles a steroid nucleus.
    # We require at least 4 rings and at least one five-membered ring.
    candidate_component = None
    for comp in components:
        if len(comp) < 4:
            continue
        # Check ring sizes in the component.
        ring_sizes = [len(rings[i]) for i in comp]
        if any(size == 5 for size in ring_sizes):
            candidate_component = comp
            break
    if candidate_component is None:
        return False, "No fused steroid nucleus (>=4 fused rings with a five-membered ring) found"
    
    # Prepare a set of all atom indices that are in any ring of the candidate fused system.
    candidate_ring_atoms = set()
    for idx in candidate_component:
        candidate_ring_atoms.update(rings[idx])
    
    # Search for a ring carbonyl group.
    # The SMARTS pattern "[R]C(=O)[R]" ensures the carbonyl carbon is in a ring.
    carbonyl_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "No ring carbonyl (3-oxo) group found"
    
    # For each match, check whether the carbonyl carbon lies in one (or more) of the candidate rings
    # and whether in that ring the carbonyl neighbor is involved in a double bond.
    # (Since the delta(1) double bond in 3-oxo steroids is located in ring A but may not connect directly 
    # to the carbonyl carbon, we loosen the requirement by looking at neighbors of the carbonyl in that ring.)
    for match in carbonyl_matches:
        # match tuple: (carbonyl carbon, oxygen, other ring carbon on the carbonyl group)
        carbonyl_idx = match[0]
        # Does the carbonyl atom belong to the candidate fused system?
        if carbonyl_idx not in candidate_ring_atoms:
            continue
        # Identify all rings (from candidate_component) that contain the carbonyl atom.
        rings_with_carbonyl = []
        for idx in candidate_component:
            if carbonyl_idx in rings[idx]:
                rings_with_carbonyl.append(rings[idx])
        # For each ring that contains the carbonyl, examine its neighbors (within the ring).
        atom_c = mol.GetAtomWithIdx(carbonyl_idx)
        for ring_atoms in rings_with_carbonyl:
            # Get neighbors in the ring (only consider non-oxygen atoms for this check).
            neighbors_in_ring = [nb.GetIdx() for nb in atom_c.GetNeighbors() if nb.GetIdx() in ring_atoms and nb.GetSymbol() == "C"]
            # For each such neighbor, check if it is engaged in a double bond with another atom also in the same ring.
            for nb_idx in neighbors_in_ring:
                nb_atom = mol.GetAtomWithIdx(nb_idx)
                for bond in nb_atom.GetBonds():
                    # Ensure the bond is double and the neighboring atom is in the same ring.
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other_idx = bond.GetOtherAtomIdx(nb_idx)
                        if other_idx in ring_atoms and other_idx != carbonyl_idx:
                            return True, "Molecule has a fused steroid nucleus with a ring 3-oxo group and a Delta(1) double bond (detected via a nearby double bond in the ring)"
    # If we did not find a suitable double bond in proximity to the carbonyl within the candidate nucleus:
    return False, "No delta(1) double bond found adjacent (within the same ring) to a 3-oxo group in the steroid nucleus"

# Example usage (uncomment to test):
# test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)C=C[C@]34C)[C@@H]1CC[C@@H]2O"  # example steroid SMILES
# result, reason = is_3_oxo_Delta_1__steroid(test_smiles)
# print(result, reason)
"""

In this implementation the criteria for a fused steroid nucleus and a nearby double bond (proxy for Delta(1)) are relaxed and expanded. This should help to avoid both false positives (non‐steroid fused systems accidentally matching the old criteria) and false negatives (steroids with extra rings or with the double bond not directly on the carbonyl carbon).