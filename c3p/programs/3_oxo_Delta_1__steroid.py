"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-Delta(1) steroid.
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
This revised implementation:
  1. Identifies fused ring components and accepts those having 3–5 rings,
     provided that at least one ring is 5-membered and at least one (or two) are 6-membered.
  2. Looks inside the candidate nucleus for a ring-carbonyl group using the SMARTS "[R]C(=O)[R]".
  3. For any ring in the nucleus that contains the carbonyl carbon, it checks for a C=C (non carbonyl)
     bond present in that same ring as a proxy for a Δ(1) double bond.
If any of these conditions are unmet, it returns False with a reason.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid.
    
    The method:
      1. Parses the molecule and obtains its rings.
      2. Builds fused ring components (rings sharing at least two atoms) and then
         accepts as candidate nuclei only those components with 3–5 fused rings that match
         a simple steroid “fingerprint”: at least one 5-membered ring and at least one six-membered ring.
      3. Within the candidate nucleus, it looks for a ring carbonyl group (using SMARTS "[R]C(=O)[R]"),
         then requires that in at least one ring (that contains that carbonyl) there is also
         at least one C=C double bond between carbon atoms (ignoring bonds that are part of the carbonyl).
    
    Args:
       smiles (str): The SMILES string of the compound.
       
    Returns:
       bool: True if the compound meets the criteria, False otherwise.
       str: A reason explaining the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of tuples of atom indices
    if not rings:
        return False, "No rings found in the molecule"
    
    n_rings = len(rings)
    # Build connectivity graph among rings: rings are fused if they share at least 2 atoms.
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(rings[i])
        for j in range(i+1, n_rings):
            set_j = set(rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Identify connected fused groups (components).
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
    
    # Define a helper to test if a fused component “looks like” a steroid nucleus.
    def is_candidate_steroid_nucleus(comp):
        # comp is a set of ring indices
        comp_rings = [rings[i] for i in comp]
        ring_sizes = [len(r) for r in comp_rings]
        count5 = sum(1 for s in ring_sizes if s == 5)
        count6 = sum(1 for s in ring_sizes if s == 6)
        # Classical steroid nucleus is 4 rings in 6-6-6-5 arrangement,
        # but we allow some variation (3 to 5 rings) to catch modified systems.
        if len(comp) < 3 or len(comp) > 5:
            return False
        if count5 < 1 or count6 < 1:
            return False
        return True

    # Look for candidate fused nucleus components.
    candidate_components = []
    for comp in components:
        if is_candidate_steroid_nucleus(comp):
            candidate_components.append(comp)
    if not candidate_components:
        return False, "No fused steroid nucleus (3–5 fused rings with at least one 5-membered and one 6-membered ring) found"
    
    # Merge candidate nucleus atoms from all candidate components.
    candidate_atoms = set()
    for comp in candidate_components:
        for idx in comp:
            candidate_atoms.update(rings[idx])
    
    # Look for a ring carbonyl group.
    # The SMARTS "[R]C(=O)[R]" forces the carbonyl carbon to be in a ring and flanked by ring atoms.
    carbonyl_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "No ring carbonyl (3-oxo) group found"
    
    # We are only interested in carbonyl hits that are inside the candidate nucleus.
    valid_carbonyl_matches = []
    for match in carbonyl_matches:
        # match is a tuple and we take the first index (the carbonyl carbon)
        if match[0] in candidate_atoms:
            valid_carbonyl_matches.append(match)
    if not valid_carbonyl_matches:
        return False, "No ring carbonyl group found within the steroid nucleus"
    
    # Identify all double bonds between carbon atoms that are not part of a C=O.
    candidate_double_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only consider C=C bonds.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                candidate_double_bonds.append((a1.GetIdx(), a2.GetIdx()))
    
    # For each valid carbonyl, find if any ring (in the candidate nucleus)
    # that contains the carbonyl atom also has a C=C double bond.
    for match in valid_carbonyl_matches:
        carbonyl_idx = match[0]
        # Identify all rings (from any candidate component) that contain this carbonyl atom.
        rings_with_carbonyl = []
        for comp in candidate_components:
            for i in comp:
                ring_atoms = rings[i]
                if carbonyl_idx in ring_atoms:
                    rings_with_carbonyl.append(ring_atoms)
        # For each such ring, check for a C=C bond.
        for ring_atoms in rings_with_carbonyl:
            l = len(ring_atoms)
            # Check bonds along the cycle (cyclically connected)
            for j in range(l):
                k = (j + 1) % l
                a1_idx = ring_atoms[j]
                a2_idx = ring_atoms[k]
                # We only consider bonds between carbons.
                bond = mol.GetBondBetweenAtoms(a1_idx, a2_idx)
                if bond is None:
                    continue
                # If this double bond is already in candidate_double_bonds, then we have a C=C.
                if (a1_idx, a2_idx) in candidate_double_bonds or (a2_idx, a1_idx) in candidate_double_bonds:
                    # We found a double bond in the same ring as a carbonyl.
                    return True, ("Molecule has a fused steroid nucleus with a ring carbonyl group and a "
                                  "C=C double bond (proxy for the Δ(1) bond) in a ring containing the carbonyl.")
    
    return False, ("No C=C double bond (other than the carbonyl’s C=O) was found in any ring containing "
                   "a ring carbonyl within the candidate steroid nucleus.")


# Example usage:
if __name__ == '__main__':
    # You can test using one of the examples.
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)C=C[C@]34C)[C@@H]1CC[C@@H]2O"  # example steroid SMILES
    result, reason = is_3_oxo_Delta_1__steroid(test_smiles)
    print(result, reason)