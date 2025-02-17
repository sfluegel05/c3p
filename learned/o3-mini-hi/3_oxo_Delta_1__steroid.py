"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-Delta(1) steroid.
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
Revised implementation based on previous feedback.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid.
    The method proceeds as follows:
      1. Parses the molecule and removes explicit hydrogens.
      2. Obtains ring information and builds fused-ring components (with rings considered fused if they share ≥2 atoms).
      3. Filters candidate steroid nuclei – they must have between 3 and 5 fused rings,
         must contain at least two 6-membered rings (as most classical steroid nuclei have three 6-membered rings),
         and the fused system should include roughly 15–22 carbons.
      4. Searches for a ring carbonyl group using the SMARTS "[R]C(=O)[R]" within the candidate nucleus.
         (This finds a carbonyl carbon that is in a ring with ring–neighbors on both sides.)
      5. For any candidate ring (within the candidate nucleus) that contains the carbonyl,
         it looks for a non–carbonyl double bond (C=C between two carbons, not involving oxygen)
         in the same ring. This double bond is used as a proxy for the Δ(1) bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria for a 3-oxo-Delta(1) steroid, False otherwise.
        str: A reason explaining the classification.
    """
    # 1. Parse molecule and remove explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.RemoveHs(mol)
    
    # 2. Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if not rings:
        return False, "No rings found in the molecule"
    n_rings = len(rings)
    
    # 3. Build connectivity graph among rings: consider rings fused if they share ≥ 2 atoms.
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(rings[i])
        for j in range(i + 1, n_rings):
            set_j = set(rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Identify fused ring components (connected components of rings).
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
    
    # Helper: determine whether a fused component “looks like” a steroid nucleus.
    def is_candidate_steroid_nucleus(comp):
        comp_rings = [rings[i] for i in comp]
        ring_sizes = [len(r) for r in comp_rings]
        # Require between 3 and 5 rings.
        if len(comp) < 3 or len(comp) > 5:
            return False
        # Require at least two 6-membered rings.
        if sum(1 for s in ring_sizes if s == 6) < 2:
            return False
        # Collect unique atoms from the fused system and count carbons.
        comp_atoms = set()
        for idx in comp:
            comp_atoms.update(rings[idx])
        carbon_count = sum(1 for a in comp_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 6)
        # Relaxed carbon count criteria (the classical steroid nucleus usually has 17 carbons,
        # so we check for roughly between 15 and 22).
        if carbon_count < 15 or carbon_count > 22:
            return False
        return True

    candidate_components = [comp for comp in components if is_candidate_steroid_nucleus(comp)]
    if not candidate_components:
        return False, ("No fused steroid nucleus (3–5 fused rings with at least two 6-membered rings and "
                       "a total of roughly 15–22 carbons) found")
    
    # Merge candidate nucleus atoms from all candidate components.
    candidate_atoms = set()
    for comp in candidate_components:
        for idx in comp:
            candidate_atoms.update(rings[idx])
    
    # 4. Look for a ring carbonyl group (3-oxo) within the candidate nucleus.
    carbonyl_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "No ring carbonyl (3-oxo) group found"
    
    # Accept only those carbonyl matches that lie within the candidate nucleus.
    valid_carbonyl_matches = []
    for match in carbonyl_matches:
        # match is a tuple; first atom is taken to be the carbonyl carbon.
        if match[0] in candidate_atoms:
            valid_carbonyl_matches.append(match)
    if not valid_carbonyl_matches:
        return False, "No ring carbonyl group found within the candidate steroid nucleus"
    
    # 5. Identify all bonds that are C=C (double bonds between two carbons).
    double_bonds = set()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bonds.add(frozenset((a1.GetIdx(), a2.GetIdx())))
    
    # For each valid carbonyl, check rings in the candidate nucleus that include its carbonyl carbon.
    for match in valid_carbonyl_matches:
        carbonyl_idx = match[0]
        # Identify all candidate nucleus rings that contain the carbonyl carbon.
        rings_with_carbonyl = []
        for comp in candidate_components:
            for i in comp:
                ring_atoms = rings[i]
                if carbonyl_idx in ring_atoms:
                    rings_with_carbonyl.append(ring_atoms)
        # In each such ring, scan each edge (bond) along the ring.
        for ring_atoms in rings_with_carbonyl:
            n = len(ring_atoms)
            for j in range(n):
                k = (j + 1) % n
                a1_idx = ring_atoms[j]
                a2_idx = ring_atoms[k]
                # Fetch the bond between these two atoms.
                bond = mol.GetBondBetweenAtoms(a1_idx, a2_idx)
                if bond is None:
                    continue
                # Skip bond if any end is oxygen (to avoid the C=O of the carbonyl).
                a1 = mol.GetAtomWithIdx(a1_idx)
                a2 = mol.GetAtomWithIdx(a2_idx)
                if a1.GetAtomicNum() == 8 or a2.GetAtomicNum() == 8:
                    continue
                # Check if this bond is a C=C double bond.
                if frozenset((a1_idx, a2_idx)) in double_bonds:
                    return True, ("Molecule has a fused steroid nucleus (3–5 fused rings with a suitable carbon count) "
                                  "with a ring carbonyl and an additional non-carbonyl C=C double bond (proxy for Δ(1)).")
    
    return False, ("No eligible C=C double bond (other than the carbonyl’s C=O) was found in any ring containing "
                   "a ring carbonyl within the candidate steroid nucleus.")


# Example usage:
if __name__ == '__main__':
    # Try with a known 3-oxo-Delta(1) steroid: helvolic acid methyl ester.
    test_smiles = "C=1[C@@]2([C@@]3(CC[C@@]/4([C@@]([C@]3(C([C@H]([C@]2([C@@H](C(C1)=O)C)[H])OC(=O)C)=O)C)(C[C@@H](\\C4=C(\\CCC=C(C)C)/C(=O)OC)OC(=O)C)C)[H])[H])C"
    result, reason = is_3_oxo_Delta_1__steroid(test_smiles)
    print(result, reason)