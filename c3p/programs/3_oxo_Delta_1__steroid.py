"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-Delta(1) steroid.
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
This implementation:
  1. Removes explicit hydrogens and obtains ring information.
  2. Builds fused ring components (rings are fused if they share ≥2 atoms),
     then accepts candidate steroid nuclei only if they have 3–5 rings,
     including at least one 5-membered and one 6-membered ring, and if the
     candidate nucleus (i.e. the union of atoms in the fused system) contains
     about 15–21 carbons (a rough estimate for the classical steroid nucleus).
  3. Searches for a ring carbonyl group within the candidate nucleus (via the SMARTS "[R]C(=O)[R]")
  4. Within any ring that contains the carbonyl atom, it requires that there is
     at least one additional C=C double bond (with both atoms being carbons) which is
     not the carbonyl’s C=O.
If all these criteria are met, the molecule is classified as a 3-oxo-Delta(1) steroid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid.
    
    The method:
      1. Parses the molecule and removes explicit hydrogens.
      2. Obtains ring information and builds fused-ring components (shared edge = ≥2 atoms).
      3. Filters candidate steroid nuclei – they must have 3–5 fused rings,
         include at least one 5-membered and one 6-membered ring, and the union of atoms
         in the nucleus should include roughly 15–21 carbons.
      4. Locates a ring carbonyl group (using SMARTS "[R]C(=O)[R]") within the candidate nucleus.
      5. For each ring (in the candidate nucleus containing the carbonyl),
         it verifies that there is at least one double bond (C=C) not belonging to the carbonyl.
    
    Args:
        smiles (str): SMILES string of the compound.
    
    Returns:
        bool: True if the compound meets the criteria, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse molecule and remove explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.RemoveHs(mol)
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if not rings:
        return False, "No rings found in the molecule"
    
    n_rings = len(rings)
    
    # Build connectivity graph among rings: consider two rings fused if they share >= 2 atoms.
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(rings[i])
        for j in range(i+1, n_rings):
            set_j = set(rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Identify connected fused groups (components) of rings.
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
    
    # Helper function: tests whether a fused component “looks like” a steroid nucleus.
    def is_candidate_steroid_nucleus(comp):
        comp_rings = [rings[i] for i in comp]
        ring_sizes = [len(r) for r in comp_rings]
        count5 = sum(1 for s in ring_sizes if s == 5)
        count6 = sum(1 for s in ring_sizes if s == 6)
        # Require between 3 and 5 rings.
        if len(comp) < 3 or len(comp) > 5:
            return False
        # Must have at least one 5-membered and one 6-membered ring.
        if count5 < 1 or count6 < 1:
            return False
        # Also check that the total number of carbons in the fused nucleus is in a reasonable range.
        comp_atoms = set()
        for idx in comp:
            comp_atoms.update(rings[idx])
        carbon_count = sum(1 for a in comp_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 6)
        if carbon_count < 15 or carbon_count > 21:
            return False
        return True

    candidate_components = []
    for comp in components:
        if is_candidate_steroid_nucleus(comp):
            candidate_components.append(comp)
    if not candidate_components:
        return False, ("No fused steroid nucleus (3–5 fused rings with at least one 5-membered "
                       "and one 6-membered ring and appropriate carbon count) found")
    
    # Merge candidate nucleus atoms from all candidate components.
    candidate_atoms = set()
    for comp in candidate_components:
        for idx in comp:
            candidate_atoms.update(rings[idx])
    
    # Look for a ring carbonyl group (3-oxo) within the candidate nucleus.
    # "[R]C(=O)[R]" forces the carbonyl carbon to be in a ring and flanked by ring atoms.
    carbonyl_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "No ring carbonyl (3-oxo) group found"
    
    # Keep only matches where the carbonyl (first atom in the match) is in the candidate nucleus.
    valid_carbonyl_matches = []
    for match in carbonyl_matches:
        if match[0] in candidate_atoms:
            valid_carbonyl_matches.append(match)
    if not valid_carbonyl_matches:
        return False, "No ring carbonyl group found within the candidate steroid nucleus"
    
    # Identify all bonds that are C=C (double bonds between two carbons)
    double_bonds = set()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bonds.add(frozenset((a1.GetIdx(), a2.GetIdx())))
    
    # For each valid carbonyl, check the rings (from candidate components)
    # that contain the carbonyl atom for a non-carbonyl C=C double bond.
    for match in valid_carbonyl_matches:
        carbonyl_idx = match[0]
        # Find rings (from any candidate candidate_component) that contain the carbonyl.
        rings_with_carbonyl = []
        for comp in candidate_components:
            for i in comp:
                ring_atoms = rings[i]
                if carbonyl_idx in ring_atoms:
                    rings_with_carbonyl.append(ring_atoms)
        # For each such ring, iterate along its (cyclic) edge bonds.
        for ring_atoms in rings_with_carbonyl:
            n = len(ring_atoms)
            for j in range(n):
                k = (j + 1) % n
                a1_idx = ring_atoms[j]
                a2_idx = ring_atoms[k]
                # Get the bond between these two atoms.
                bond = mol.GetBondBetweenAtoms(a1_idx, a2_idx)
                if bond is None:
                    continue
                # Skip the bond if it involves an oxygen (i.e. is likely part of the carbonyl)
                a1 = mol.GetAtomWithIdx(a1_idx)
                a2 = mol.GetAtomWithIdx(a2_idx)
                if a1.GetAtomicNum() == 8 or a2.GetAtomicNum() == 8:
                    continue
                # Check if this bond is a C=C double bond.
                if frozenset((a1_idx, a2_idx)) in double_bonds:
                    return True, ("Molecule has a fused steroid nucleus with a ring carbonyl group and a "
                                  "C=C double bond (proxy for the Δ(1) bond) in a ring containing the carbonyl.")
    
    return False, ("No C=C double bond (other than the carbonyl’s C=O) was found in any ring containing "
                   "a ring carbonyl within the candidate steroid nucleus.")

# Example usage:
if __name__ == '__main__':
    # Test with one example (e.g. helvolic acid methyl ester)
    test_smiles = "C=1[C@@]2([C@@]3(CC[C@@]/4([C@@]([C@]3(C([C@H]([C@]2([C@@H](C(C1)=O)C)[H])OC(=O)C)=O)C)(C[C@@H](\\C4=C(\\CCC=C(C)C)/C(=O)OC)OC(=O)C)C)[H])[H])C"
    result, reason = is_3_oxo_Delta_1__steroid(test_smiles)
    print(result, reason)