"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-Delta(1) steroid.
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
This implementation:
  1. Identifies a fused steroid nucleus by grouping rings that share at least two atoms.
     A candidate nucleus requires ≥4 rings (often 3 six-membered and 1 five-membered ring).
  2. Looks for a ring carbonyl group using the SMARTS pattern "[R]C(=O)[R]".
  3. For each candidate ring (i.e. ring in the fused component that contains the carbonyl),
     it scans the ordered ring (from RDKit’s GetRingInfo) for a carbon–carbon double bond,
     which stands in as a proxy for a Delta(1) double bond.
If any of these conditions are unmet, it returns False with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid.
    
    It first attempts to detect a steroid nucleus as a fused ring system (≥4 rings with at least one five-membered ring).
    Then it searches for a ring carbonyl ([R]C(=O)[R]) within that nucleus.
    Finally, it checks that in at least one ring containing the carbonyl, there exists a double bond
    (between two carbons) that is not part of the carbonyl group—serving as a proxy for the Δ(1) (1,2) double bond.

    Args:
        smiles (str): The SMILES string of the compound.
        
    Returns:
        bool: True if the molecule meets the criteria for a 3-oxo-Delta(1) steroid, False otherwise.
        str: A reason explaining the result.
    """
    
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # Build a connectivity graph between rings.
    n_rings = len(rings)
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(rings[i])
        for j in range(i+1, n_rings):
            set_j = set(rings[j])
            if len(set_i.intersection(set_j)) >= 2:  # fused if share 2 or more atoms
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Identify connected components of rings.
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
    
    # Look for a candidate fused steroid nucleus: one component with ≥4 rings and at least one 5-membered ring.
    candidate_component = None
    for comp in components:
        if len(comp) < 4:
            continue
        ring_sizes = [len(rings[i]) for i in comp]
        if any(size == 5 for size in ring_sizes):
            candidate_component = comp
            break
    if candidate_component is None:
        return False, "No fused steroid nucleus (≥4 fused rings with at least one 5-membered ring) found"
    
    # Build a set comprising all atom indices in the candidate fused nucleus.
    candidate_atoms = set()
    for idx in candidate_component:
        candidate_atoms.update(rings[idx])
    
    # First, search for a ring carbonyl group.
    # The SMARTS "[R]C(=O)[R]" ensures the carbonyl carbon is in a ring and flanked by ring atoms.
    carbonyl_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "No ring carbonyl (3-oxo) group found"
    
    # Next, search for a C=C double bond between two ring carbons (ignoring the carbonyl C=O).
    # We require that at least one of the rings (from candidate_component) that contains a carbonyl
    # also has a C=C double bond.
    # Get all bonds in the molecule that are double and between carbon atoms.
    candidate_double_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                candidate_double_bonds.append((a1.GetIdx(), a2.GetIdx()))
    
    # For each carbonyl match, check if it is part of a candidate ring that also has a double bond.
    for match in carbonyl_matches:
        carbonyl_idx = match[0]  # the carbonyl carbon index
        # Only consider if the carbonyl is in the candidate nucleus.
        if carbonyl_idx not in candidate_atoms:
            continue
        
        # Identify all rings (by index from candidate_component) that contain the carbonyl atom.
        rings_with_carbonyl = [rings[i] for i in candidate_component if carbonyl_idx in rings[i]]
        # For each such ring, check if one of its bonds (in order) is a double bond between two carbons.
        for ring in rings_with_carbonyl:
            # 'ring' is a tuple that represents an ordered cycle.
            l = len(ring)
            # Check each bond in the cycle 
            for i in range(l):
                j = (i + 1) % l  # next atom in the cycle (cyclic)
                a1_idx = ring[i]
                a2_idx = ring[j]
                # Skip if this bond belongs to the carbonyl pattern (i.e. from the carbonyl C to an O).
                # Our candidate_double_bonds list only has C-C double bonds so we are safe.
                if (a1_idx, a2_idx) in candidate_double_bonds or (a2_idx, a1_idx) in candidate_double_bonds:
                    # Found a qualifying double bond.
                    return True, ("Molecule has a fused steroid nucleus with a ring carbonyl group and a C=C double bond "
                                  "in the same ring (proxy for the Delta(1) double bond).")
    
    return False, ("No C=C double bond (other than the carbonyl's C=O) was found in any ring containing a ring carbonyl "
                   "within the fused steroid nucleus.")
   
# Example usage (uncomment to test):
# test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)C=C[C@]34C)[C@@H]1CC[C@@H]2O"  # example steroid SMILES
# result, reason = is_3_oxo_Delta_1__steroid(test_smiles)
# print(result, reason)