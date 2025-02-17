"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid
Definition: Any terpenoid derived from a triterpene (typically built from a C30 skeleton that may be rearranged or modified).
Heuristic improvements include:
  - Molecular weight > 400 Da
  - Overall high fraction of sp3 carbons (> 0.5)
  - Limited heteroatom (especially oxygen) content relative to carbons (to avoid highly decorated sugar‐moieties)
  - Extraction of rings: if at least 4 predominantly non‐aromatic (aromatic atom fraction < 0.3) rings exist, use these to identify a fused core.
    Else, use all rings.
  - Build connected (fused) clusters and choose the largest union as the candidate core.
  - Expand the core by iteratively adding neighboring carbon atoms that attach to at least two atoms of the core.
  - Finally require that the expanded core have a carbon count between 15 and 35, and that few of these carbons are aromatic.
  
The additional oxygen-to-carbon ratio check helps reject glycosylated or sugar‐rich compounds.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string
    using improved heuristics.

    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if classified as a triterpenoid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Overall molecular weight must be sufficiently high for a triterpenoid.
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a triterpenoid"
    
    # Overall fraction of sp3 carbons should be high (triterpenoids are usually saturated).
    frac_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    if frac_sp3 < 0.5:
        return False, f"Low fraction of sp3 carbons ({frac_sp3:.2f}); triterpenoids are usually highly saturated"
    
    # Check heteroatom content vs. carbon: too many non-carbon heavy atoms disfavors a triterpenoid core.
    atoms = mol.GetAtoms()
    n_carbons = sum(1 for a in atoms if a.GetAtomicNum() == 6)
    # count oxygen atoms separately (they often come from sugars)
    n_oxygens = sum(1 for a in atoms if a.GetAtomicNum() == 8)
    n_heavy = sum(1 for a in atoms if a.GetAtomicNum() > 1)
    if n_heavy > 0 and ((n_heavy - n_carbons) / n_heavy) > 0.4:
        return False, f"High fraction of heteroatoms ({(n_heavy - n_carbons) / n_heavy:.2f}); not typical for triterpenoids"
    # Extra check: if oxygen/carbon ratio is too high, suspect sugar moieties.
    if n_carbons > 0 and (n_oxygens / n_carbons) > 0.35:
        return False, f"High oxygen-to-carbon ratio ({n_oxygens}/{n_carbons}); likely contains sugar moieties"

    # Get ring information.
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings detected in the molecule."
    
    # Identify rings that are predominantly non‐aromatic.
    non_aromatic_rings = []
    for ring in all_rings:
        aromatic_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
        if (aromatic_count / len(ring)) < 0.3:
            non_aromatic_rings.append(ring)
    
    # Try to use non‐aromatic rings if there are at least four; otherwise fall back on all rings.
    if len(non_aromatic_rings) >= 4:
        candidate_rings = non_aromatic_rings
        used_rings_info = f"{len(candidate_rings)} non‐aromatic rings"
    else:
        candidate_rings = all_rings
        used_rings_info = f"all rings (total {len(candidate_rings)})"
    
    # Build connectivity among candidate rings: two rings are connected if they share an atom.
    n_rings = len(candidate_rings)
    adjacency = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if set(candidate_rings[i]).intersection(candidate_rings[j]):
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Find connected components (clusters) among the candidate rings.
    visited = set()
    components = []
    for i in range(n_rings):
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
    
    # Choose the fused cluster that gives the maximum union of atoms.
    best_union = set()
    best_component = None
    for comp in components:
        union_atoms = set()
        for ring_index in comp:
            union_atoms.update(candidate_rings[ring_index])
        if len(union_atoms) > len(best_union):
            best_union = union_atoms
            best_component = comp

    # If our chosen cluster has at least four rings, use it.
    if best_component is not None and len(best_component) >= 4:
        core_atoms = set(best_union)
        cluster_info = f"Fused cluster with {len(best_component)} rings"
    else:
        # Otherwise, fall back on the union of candidate rings.
        core_atoms = set()
        for ring in candidate_rings:
            core_atoms.update(ring)
        cluster_info = f"Union of candidate rings ({len(candidate_rings)} rings)"
    
    # Expand the core by adding any carbon atom not in the core that has at least 2 neighbors already in the core.
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
    
    # Compute the fraction of aromatic carbons in the expanded core.
    aromatic_core = sum(1 for idx in core_carbon_idxs if mol.GetAtomWithIdx(idx).GetIsAromatic())
    frac_aromatic_core = aromatic_core / n_core_carbons if core_carbon_idxs else 0.0
    
    # For many true triterpenes, the core (expanded fused ring system) has between 15 and 35 carbons.
    if not (15 <= n_core_carbons <= 35):
        return False, (f"Expanded core has {n_core_carbons} carbons (expected between 15 and 35) "
                       f"using {cluster_info} from {used_rings_info}")
    
    # A high aromatic fraction in the core is not expected.
    if frac_aromatic_core >= 0.3:
        return False, f"Aromatic fraction in core is {frac_aromatic_core:.2f} (expected below 0.3)"
    
    return True, (f"{cluster_info} expanded to a core with {n_core_carbons} carbons, "
                  f"aromatic fraction: {frac_aromatic_core:.2f}, overall sp3 fraction: {frac_sp3:.2f}, "
                  "is consistent with a triterpenoid.")

# For demonstration or testing:
if __name__ == "__main__":
    # Example: Deacetylnomilin (a known triterpenoid)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)