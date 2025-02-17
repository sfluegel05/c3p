"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5β-steroid
Definition: Any 3-oxo steroid that has beta- configuration at position 5.
This implementation uses the following heuristics:
  • It requires the presence of a fused, non‐aromatic ring system formed by at least 4 rings,
    whose union contains at least 17 carbon atoms (the typical steroid nucleus).
  • It requires a ring-bound ketone group (the 3-oxo moiety) using the SMARTS "[R]C(=O)[R]".
  • It requires that the isomeric SMILES string includes at least one '@@' as a proxy for 5β configuration.
Note: A completely correct detection would require explicit identification of the steroid numbering system.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5β-steroid.
    
    Heuristic criteria:
      1. The molecule must have a fused, non-aromatic ring system corresponding to a steroid nucleus.
         In this implementation, we:
           - Find all non-aromatic rings.
           - Build a graph of rings that share at least 2 atoms (fused rings).
           - Require that the largest connected (fused) ring system is composed of at least 4 rings,
             and that the union of atoms in these rings includes at least 17 carbon atoms.
      2. The molecule must contain at least one ring-bound ketone group (the 3-oxo moiety), detected via
         the SMARTS "[R]C(=O)[R]".
      3. The molecule’s isomeric SMILES must include at least one chiral specification with '@@',
         taken as a proxy for the 5β configuration.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a 3-oxo-5β-steroid, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that stereochemistry is assigned.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # --- Criterion 1: Fused steroid nucleus detection ---
    ring_info = mol.GetRingInfo()
    ring_atom_sets = []
    # Get only non-aromatic rings
    for ring in ring_info.AtomRings():
        if all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            ring_atom_sets.append(set(ring))
    if not ring_atom_sets:
        return False, "No non-aromatic rings found; no steroid nucleus"
    
    # Build a connectivity graph: nodes are rings and an edge exists if two rings share at least 2 atoms.
    n = len(ring_atom_sets)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(ring_atom_sets[i].intersection(ring_atom_sets[j])) >= 2:
                graph[i].add(j)
                graph[j].add(i)
    
    # Find connected components (fused ring systems)
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
    
    # Evaluate each fused component for the steroid nucleus: at least 4 rings and at least 17 carbons.
    steroid_nucleus_found = False
    for comp in fused_components:
        if len(comp) >= 4:
            # Get union of atoms in these rings
            atom_union = set()
            for idx in comp:
                atom_union = atom_union.union(ring_atom_sets[idx])
            # Count carbon atoms in the union
            carbon_count = sum(1 for atom_idx in atom_union if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
            if carbon_count >= 17:
                steroid_nucleus_found = True
                break
    if not steroid_nucleus_found:
        return False, "No fused ring system with at least 4 rings and 17 carbons (steroid nucleus) found"
    
    # --- Criterion 2: 3-oxo (ring-bound ketone) ---
    # SMARTS for a ketone with both neighbors in a ring: "[R]C(=O)[R]"
    ketone_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ring-bound ketone (3-oxo) group found"
    
    # --- Criterion 3: Proxy for 5β configuration ---
    iso_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    if "@@" not in iso_smi:
        return False, "No chiral center with '@@' (indicative of beta configuration) detected"
    
    return True, "Molecule has a fused non-aromatic ring system (>=4 rings, >=17 carbons), a ring-bound ketone (3-oxo), and a '@@' chiral center indicative of 5β configuration"

# Example usage (uncomment to test):
#if __name__ == "__main__":
#    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2C(=O)CO"  # 5β-dihydrodeoxycorticosterone
#    result, reason = is_3_oxo_5beta_steroid(test_smiles)
#    print(result, reason)