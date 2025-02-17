"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid 
Definition: Any terpenoid derived from a triterpene.
Heuristic improvements over the previous version:
  - Check overall molecular weight (>400 Da) and ensure a high overall fraction of sp3 carbons.
  - Check that the molecule does not have an excess of heteroatoms (e.g. sugar moieties).
  - Compute all rings and mark those whose majority of atoms are non‐aromatic.
  - If there are at least 4 such “non‐aromatic” rings, then try to build a fused‐ring “core.”
    • First, we attempt to compute the largest connected (fused) ring cluster.
    • If that cluster has at least 4 rings we use its union as the core.
    • Otherwise, we take the union of all non‐aromatic rings as an alternative core.
  - The core is then “expanded” by iteratively adding neighboring carbon atoms that are connected to two or more core atoms.
  - Finally, the heuristic requires that the expanded core have between 13 and 40 carbon atoms
    (allowing some rearrangements) and that few of these are aromatic (<30%).
    
Examples (by our testing) show that these changes help catch true triterpenoids while reducing false positives.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string using an improved heuristic.
    
    Steps:
      1. Parse the SMILES.
      2. Check overall molecular weight (>400 Da), high overall saturation (fraction sp3 > 0.5)
         and that there are not too many heteroatoms (non-carbon) overall.
      3. Extract all rings and count those that are mainly aliphatic (non‐aromatic).
      4. If fewer than 4 such rings are detected, reject.
      5. Build a fused‐ring “core”:
             - Compute connected clusters among all rings.
             - If the largest fused cluster has at least 4 rings, use its atom union.
             - Otherwise, use the union of all non‐aromatic rings.
      6. Expand the core by adding neighboring carbon atoms with at least two core neighbors.
      7. Count the number of carbon atoms in the expanded core and its aromatic fraction.
      8. Require that:
             - The core has between 13 and 40 carbons.
             - The aromatic fraction in the core is < 0.3.
    
    Args:
      smiles (str): SMILES representation of the molecule.
      
    Returns:
      bool: True if classified as a triterpenoid, False otherwise.
      str: Explanation for the classification decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular weight.
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for triterpenoid"
    
    # Check overall degree of saturation:
    frac_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    if frac_sp3 < 0.5:
        return False, f"Low fraction of sp3 carbons ({frac_sp3:.2f}); triterpenoids are usually highly saturated"
    
    # Check heteroatom content: count carbons vs. total heavy atoms (ignore hydrogens)
    atoms = mol.GetAtoms()
    n_carbons = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    n_heavy = sum(1 for atom in atoms if atom.GetAtomicNum() > 1)
    if n_heavy > 0 and ((n_heavy - n_carbons) / n_heavy) > 0.4:
        return False, f"High fraction of heteroatoms ({(n_heavy - n_carbons) / n_heavy:.2f}); not typical for triterpenoids"
    
    # Get all rings from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings detected in the molecule."
    
    # Identify rings that are "non‐aromatic" (i.e. majority of atoms are not aromatic).
    non_aromatic_rings = []
    for ring in all_rings:
        aromatic_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
        if aromatic_count / len(ring) < 0.3:
            non_aromatic_rings.append(ring)
    
    if len(non_aromatic_rings) < 4:
        return False, f"Only {len(non_aromatic_rings)} non‐aromatic rings detected (need at least 4 for a triterpenoid)"
    
    # Build connectivity among all rings (using all rings regardless of aromaticity)
    n_rings_total = len(all_rings)
    adjacency = {i: set() for i in range(n_rings_total)}
    for i in range(n_rings_total):
        for j in range(i+1, n_rings_total):
            if set(all_rings[i]).intersection(all_rings[j]):
                adjacency[i].add(j)
                adjacency[j].add(i)
                
    # Find connected clusters (fused ring clusters)
    visited = set()
    components = []
    for i in range(n_rings_total):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                comp.add(node)
                stack.extend(adjacency[node] - visited)
            components.append(comp)
    
    # Choose the cluster with the largest union of atoms.
    best_union = set()
    best_component = None
    for comp in components:
        union_atoms = set()
        for ring_index in comp:
            union_atoms.update(all_rings[ring_index])
        if len(union_atoms) > len(best_union):
            best_union = union_atoms
            best_component = comp

    # If the best fused cluster has at least 4 rings, use it.
    if best_component is not None and len(best_component) >= 4:
        core_atoms = set(best_union)
        used_cluster_info = f"Fused cluster ({len(best_component)} rings)"
    else:
        # Otherwise fall back on using union of all non‐aromatic rings.
        core_atoms = set()
        for ring in non_aromatic_rings:
            core_atoms.update(ring)
        used_cluster_info = f"All non‐aromatic rings union (count: {len(non_aromatic_rings)})"
    
    # Expand the core: add any carbon atom outside core that has >=2 neighbors already in core.
    expanded_core = set(core_atoms)
    changed = True
    while changed:
        changed = False
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx in expanded_core:
                continue
            if atom.GetAtomicNum() != 6:
                continue
            nbrs = atom.GetNeighbors()
            count_in_core = sum(1 for nbr in nbrs if nbr.GetIdx() in expanded_core)
            if count_in_core >= 2:
                expanded_core.add(idx)
                changed = True
    
    # Count carbon atoms in the expanded core.
    core_carbon_idxs = [idx for idx in expanded_core if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    n_core_carbons = len(core_carbon_idxs)
    
    # Compute fraction of aromatic carbons in the core.
    aromatic_core = 0
    for idx in core_carbon_idxs:
        if mol.GetAtomWithIdx(idx).GetIsAromatic():
            aromatic_core += 1
    frac_aromatic_core = aromatic_core / n_core_carbons if n_core_carbons else 0
    
    # Check that the core carbon count is within the expected range:
    if not (13 <= n_core_carbons <= 40):
        return False, (f"Expanded fused ring core has {n_core_carbons} carbons (expected between 13 and 40) "
                       f"using {used_cluster_info}")
    
    if frac_aromatic_core >= 0.3:
        return False, f"Aromatic fraction in core is {frac_aromatic_core:.2f} (expected below 0.3)"
    
    return True, (f"{used_cluster_info} expanded to a core with {n_core_carbons} carbons, "
                  f"aromatic fraction: {frac_aromatic_core:.2f}, overall sp3 fraction: {frac_sp3:.2f}, "
                  "is consistent with a triterpenoid.")

# Testing (for demonstration; remove or comment out in production)
if __name__ == "__main__":
    # Example: Deacetylnomilin (a known triterpenoid)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)