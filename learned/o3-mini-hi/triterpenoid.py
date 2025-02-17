"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid
Definition: Any terpenoid derived from a triterpene (typically built from a C30 skeleton that may be rearranged or modified).
Heuristic improvements include:
  - Removal of sugar (glycoside) moieties prior to analysis.
  - Molecular weight > 400 Da (after sugar removal) 
  - Overall high fraction of sp3 carbons (> 0.5)
  - Limited heteroatom content (especially oxygen) relative to carbons, to avoid sugar‐rich compounds.
  - Extraction of rings: if at least 4 predominantly non‐aromatic (aromatic fraction < 0.3) rings exist, use these to identify a fused core.
    Else, use all rings.
  - Build connected (fused) clusters and choose the largest union as the candidate core.
  - Expand the core by iteratively adding neighboring carbon atoms that attach to at least two atoms of the core.
  - Finally require that the expanded core have a carbon count between 15 and 35 and that few of these carbons are aromatic.
  
Sugar removal: we “peel off” rings that are 5 or 6 members and that have a high oxygen content (>=50% O atoms). 
This helps avoid classifying a true triterpenoid as “non‐triterpenoid” because it is decorated with sugars, and 
also helps avoid false positives on compounds that are dominated by sugar features.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string using improved heuristics.
    First attempts to remove sugar moieties (based on detecting oxygen-rich 5–6 membered rings)
    and then applies several tests on the remaining core structure.

    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if classified as a triterpenoid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- STEP 1: Remove likely sugar rings to avoid high O/C artifacts -----
    # We inspect all rings (using the ring info) and mark 5- or 6-membered rings with >=50% oxygen atoms.
    ring_info = mol.GetRingInfo()
    sugar_atoms = set()
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            nO = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if nO / len(ring) >= 0.5:
                sugar_atoms.update(ring)
    # If any sugar candidates were found, remove them.
    if sugar_atoms:
        # Create an editable copy and remove atoms in descending order (to preserve atom indices)
        rw_mol = Chem.RWMol(mol)
        for idx in sorted(sugar_atoms, reverse=True):
            # Removal of an atom will remove attached bonds.
            rw_mol.RemoveAtom(idx)
        mol_no_sugar = rw_mol.GetMol()
        # Sometimes removal leads to disconnected fragments; choose the largest fragment.
        frags = Chem.GetMolFrags(mol_no_sugar, asMols=True)
        if not frags:
            return False, "No fragments remain after sugar removal"
        mol = max(frags, key=lambda m: m.GetNumAtoms())
        sugar_info = " (sugar moieties removed)"
    else:
        sugar_info = ""
    
    # ----- STEP 2: Basic molecular properties -----
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a triterpenoid{sugar_info}"
    
    frac_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    if frac_sp3 < 0.5:
        return False, f"Low fraction of sp3 carbons ({frac_sp3:.2f}); triterpenoids are usually highly saturated{sugar_info}"
    
    atoms = mol.GetAtoms()
    n_carbons = sum(1 for a in atoms if a.GetAtomicNum() == 6)
    n_oxygens = sum(1 for a in atoms if a.GetAtomicNum() == 8)
    n_heavy = sum(1 for a in atoms if a.GetAtomicNum() > 1)
    if n_heavy > 0 and ((n_heavy - n_carbons) / n_heavy) > 0.4:
        return False, (f"High fraction of heteroatoms ({(n_heavy - n_carbons) / n_heavy:.2f}); " 
                       f"not typical for triterpenoids{sugar_info}")
    if n_carbons > 0 and (n_oxygens / n_carbons) > 0.35:
        return False, f"High oxygen-to-carbon ratio ({n_oxygens}/{n_carbons}); likely contains extra decorations{sugar_info}"
    
    # ----- STEP 3: Identifying candidate rings in the molecule -----
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings detected in the molecule" + sugar_info
    
    non_aromatic_rings = []
    for ring in all_rings:
        aromatic_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
        if (aromatic_count / len(ring)) < 0.3:
            non_aromatic_rings.append(ring)
    
    if len(non_aromatic_rings) >= 4:
        candidate_rings = non_aromatic_rings
        used_rings_info = f"{len(candidate_rings)} non‐aromatic rings"
    else:
        candidate_rings = all_rings
        used_rings_info = f"all rings (total {len(candidate_rings)})"
    
    # Build connectivity among candidate rings: two rings are connected if they share at least one atom.
    n_rings = len(candidate_rings)
    adjacency = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if set(candidate_rings[i]).intersection(candidate_rings[j]):
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Identify connected components (clusters) among candidate rings.
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
    
    # Choose the component (fused cluster) that gives the maximum union of ring atoms.
    best_union = set()
    best_component = None
    for comp in components:
        union_atoms = set()
        for ring_index in comp:
            union_atoms.update(candidate_rings[ring_index])
        if len(union_atoms) > len(best_union):
            best_union = union_atoms
            best_component = comp

    if best_component is not None and len(best_component) >= 4:
        core_atoms = set(best_union)
        cluster_info = f"Fused cluster with {len(best_component)} rings"
    else:
        core_atoms = set()
        for ring in candidate_rings:
            core_atoms.update(ring)
        cluster_info = f"Union of candidate rings ({len(candidate_rings)} rings)"
    
    # Expand the core by adding any carbon atom (not already in the core) that has at least 2 neighbors in the core.
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

    # Analyze the expanded core.
    core_carbon_idxs = [idx for idx in expanded_core if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    n_core_carbons = len(core_carbon_idxs)
    aromatic_core = sum(1 for idx in core_carbon_idxs if mol.GetAtomWithIdx(idx).GetIsAromatic())
    frac_aromatic_core = aromatic_core / n_core_carbons if core_carbon_idxs else 0.0

    if not (15 <= n_core_carbons <= 35):
        return False, (f"Expanded core has {n_core_carbons} carbons (expected between 15 and 35) "
                       f"using {cluster_info} from {used_rings_info}{sugar_info}")
    if frac_aromatic_core >= 0.3:
        return False, f"Aromatic fraction in core is {frac_aromatic_core:.2f} (expected below 0.3){sugar_info}"
    
    return True, (f"{cluster_info} expanded to a core with {n_core_carbons} carbons, aromatic fraction: {frac_aromatic_core:.2f}, "
                  f"overall sp3 fraction: {frac_sp3:.2f}, is consistent with a triterpenoid{sugar_info}")

# For demonstration or testing:
if __name__ == "__main__":
    # Example: Deacetylnomilin (a known triterpenoid)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)