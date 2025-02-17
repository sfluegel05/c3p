"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5β-steroid
Definition: Any 3-oxo steroid that has beta-configuration at position 5.
Revised Heuristics:
  • The molecule must have a fused non-aromatic ring system in which at least one subset of 4 rings 
    corresponds to a steroid nucleus (one 5-membered and three 6-membered rings) and that, taken together, 
    comprise at least 17 carbon atoms.
  • The molecule must contain at least one ring-bound ketone group (3-oxo), detected via SMARTS "[R]C(=O)[R]".
  • The molecule’s isomeric SMILES must include at least one '@@' (as a proxy for 5β configuration).
Note: This is only a heuristic detection; no attempt is made to fully number the steroid skeleton.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5β-steroid based on revised heuristics.
    
    Heuristic criteria:
      1. Must contain a fused, non-aromatic ring system in which there is at least one subset of 4
         rings having sizes [5,6,6,6] (the typical steroid nucleus) and containing at least 17 carbons.
      2. Must contain a ring-bound ketone group (3-oxo) detected by the SMARTS "[R]C(=O)[R]".
      3. Must have at least one '@@' chiral specification in its isomeric SMILES, taken as a proxy for 
         5β configuration.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a 3-oxo-5β-steroid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # --- Criterion 1: Find a fused ring system matching steroid nucleus ---
    ring_info = mol.GetRingInfo()
    # Only consider non-aromatic rings
    ring_atom_sets = []
    ring_sizes = []  # record ring size for each ring
    for ring in ring_info.AtomRings():
        if all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            ring_atom_sets.append(set(ring))
            ring_sizes.append(len(ring))
    if not ring_atom_sets:
        return False, "No non-aromatic rings found; no steroid nucleus"
    
    # Build a connectivity graph on rings: nodes=rings; an edge if rings share at least 2 atoms (fused)
    n = len(ring_atom_sets)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(ring_atom_sets[i].intersection(ring_atom_sets[j])) >= 2:
                graph[i].add(j)
                graph[j].add(i)
                
    # Identify connected components (fused ring systems)
    visited = set()
    fused_components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                curr = stack.pop()
                if curr not in visited:
                    visited.add(curr)
                    comp.add(curr)
                    stack.extend(graph[curr] - visited)
            fused_components.append(comp)
    
    steroid_nucleus_found = False
    # For each fused component, try to find a subset of 4 rings that match sizes [5,6,6,6].
    for comp in fused_components:
        if len(comp) < 4:
            continue  # too few rings to be the steroid nucleus
        # Get list of rings and ring sizes from this component.
        comp_ring_indices = list(comp)
        comp_ring_sizes = [ring_sizes[i] for i in comp_ring_indices]
        # Check all combinations of 4 rings from the component to see if one set has one 5-membered and three 6-membered rings.
        from itertools import combinations
        for subset in combinations(comp_ring_indices, 4):
            subset_sizes = [ring_sizes[i] for i in subset]
            # Sort for easier comparison
            sorted_sizes = sorted(subset_sizes)
            if sorted_sizes == [5, 6, 6, 6]:
                # Also check that the union of these rings covers >=17 carbon atoms.
                atom_union = set()
                for i in subset:
                    atom_union |= ring_atom_sets[i]
                carbon_count = sum(1 for idx in atom_union if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                if carbon_count >= 17:
                    steroid_nucleus_found = True
                    break
        if steroid_nucleus_found:
            break
    
    if not steroid_nucleus_found:
        return False, "No fused ring system with a steroid nucleus signature (one 5-membered and three 6-membered rings, ≥17 carbons) found"
    
    # --- Criterion 2: Find a ring-bound ketone
    ketone_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ring-bound ketone (3-oxo) group found"
    
    # --- Criterion 3: Check for beta configuration proxy using '@@'
    iso_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    if "@@" not in iso_smi:
        return False, "No chiral center with '@@' (indicative of beta configuration) detected"
    
    return True, ("Molecule has a fused steroid nucleus (subset of 4 rings with sizes [5,6,6,6] and ≥17 carbons), "
                  "a ring-bound ketone (3-oxo), and a '@@' chiral center indicative of 5β configuration")
    
# Example usage (uncomment for testing):
#if __name__ == "__main__":
#    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2C(=O)CO"  # 5β-dihydrodeoxycorticosterone
#    result, reason = is_3_oxo_5beta_steroid(test_smiles)
#    print(result, reason)