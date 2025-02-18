"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid (any terpenoid derived from a triterpene)
Heuristic:
  - Check that overall molecular weight > 400 Da.
  - Enforce that the molecule is largely aliphatic (e.g. high fraction of sp3 carbons globally).
  - Extract all rings, then build a fused ring cluster.
  - Expand the fused ring “core” by including neighboring carbon atoms connected to at least two core atoms.
  - Count the number of fused rings, the number of carbons in the core, and the aromatic fraction within that core.
  - Require that there are at least 4 fused rings, that the expanded core has a carbon count between 15 and 40,
    and that very few of the core carbons are aromatic (<30%).
This heuristic was designed after noting that some triterpenoids (especially those derived from a C30 skeleton)
can lose methyl groups or be rearranged, so the core might not have exactly 27–35 carbons.
Also, many false‐positive non‐terpenoids had many aromatic atoms.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string using an improved heuristic.
    
    Steps:
      1. Parse the SMILES.
      2. Check overall molecular weight (>400 Da).
      3. Check global saturation by computing fraction of sp3 carbons (should be >0.5).
      4. Extract all rings and compute connectivity among them. Identify the largest fused ring cluster.
      5. Expand the core of the fused-ring cluster by including neighboring carbon atoms that are connected
         to at least two atoms already in the core.
      6. Count the number of fused rings and the number of carbon atoms present in the expanded core.
      7. Also calculate the fraction of aromatic carbon atoms in the core.
      8. Require:
           - At least 4 fused rings.
           - The core has between 15 and 40 carbon atoms.
           - The core’s aromatic carbon fraction is below 0.3.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a triterpenoid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check that overall molecular weight is high enough.
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for triterpenoid"
    
    # Check overall degree of saturation.
    frac_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    if frac_sp3 < 0.5:
        return False, f"Low fraction of sp3 carbons ({frac_sp3:.2f}); triterpenoids are typically largely saturated"
    
    # Get ring info: list of tuples of atom indices.
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings detected in the molecule."
    
    n_rings_total = len(all_rings)
    
    # Build a connectivity among rings: adjacent if they share at least one atom.
    adjacency = {i: set() for i in range(n_rings_total)}
    for i in range(n_rings_total):
        for j in range(i + 1, n_rings_total):
            if set(all_rings[i]).intersection(all_rings[j]):
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Identify connected clusters of rings.
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
    
    # Choose the fused ring cluster whose union of atoms is largest.
    best_union = set()
    best_component = None
    for comp in components:
        union_atoms = set()
        for ring_index in comp:
            union_atoms.update(all_rings[ring_index])
        if len(union_atoms) > len(best_union):
            best_union = union_atoms
            best_component = comp
    
    if best_component is None:
        return False, "No fused ring cluster detected."
    
    n_fused_rings = len(best_component)
    
    # Expand the fused ring core by iteratively adding carbon atoms outside the union if 
    # they have at least two neighbors already in the core.
    core = set(best_union)
    added = True
    while added:
        added = False
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx in core:
                continue
            # Consider only carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Count neighbors in the core.
            nbrs = atom.GetNeighbors()
            neighbor_count = sum(1 for nbr in nbrs if nbr.GetIdx() in core)
            if neighbor_count >= 2:
                core.add(idx)
                added = True
    
    # Count number of carbon atoms in core.
    core_carbons = [idx for idx in core if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    n_core_carbons = len(core_carbons)
    
    # Also, count how many of these core carbons are aromatic.
    aromatic_core_carbons = 0
    for idx in core_carbons:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetIsAromatic():
            aromatic_core_carbons += 1
    if n_core_carbons > 0:
        frac_aromatic_core = aromatic_core_carbons / n_core_carbons
    else:
        frac_aromatic_core = 0
        
    # Apply heuristic criteria.
    if n_fused_rings < 4:
        return False, f"Fused ring cluster has only {n_fused_rings} rings (expected at least 4 for triterpenoid)"
    # Allow a wider range for modified/rearranged cores.
    if not (15 <= n_core_carbons <= 40):
        return False, f"Expanded fused ring core has {n_core_carbons} carbons (expected between 15 and 40 for triterpenoid core)"
    if frac_aromatic_core >= 0.3:
        return False, f"Aromatic fraction in core is {frac_aromatic_core:.2f} (expected below 0.3 since triterpenoids are mostly aliphatic)"
    
    return True, (f"Fused ring core with {n_fused_rings} rings and {n_core_carbons} carbons, "
                  f"aromatic fraction in core: {frac_aromatic_core:.2f}, overall sp3 fraction: {frac_sp3:.2f}, "
                  "is consistent with a triterpenoid.")

# For testing purposes (this section may be removed or commented out in production)
if __name__ == "__main__":
    # Example: Deacetylnomilin (a triterpenoid) SMILES string
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)