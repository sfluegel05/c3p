"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5β-steroid
Definition: Any 3-oxo steroid that has beta-configuration at position 5.
Revised Heuristics:
  • Identify fused, non-aromatic ring systems and look for a candidate steroid nucleus comprised
    of exactly 4 rings with sizes [5, 6, 6, 6] (sorted order) where the union of rings contains ≥17 C atoms.
  • Require that at least one ring-bound ketone (detected via SMARTS "[R]C(=O)[R]") is located
    within the candidate nucleus.
  • Require that the isomeric SMILES contains at least one '@@', taken as a proxy for 5β configuration.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from itertools import combinations

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5β-steroid based on revised heuristics.
    
    Heuristic criteria:
      1. Must contain a fused, non-aromatic ring system that includes a subset of 4 rings whose sizes,
         when sorted, equal [5,6,6,6] and the union of these rings covers ≥17 carbon atoms.
      2. Must contain a ring-bound ketone group (3-oxo) that is part of the candidate nucleus.
      3. Must have at least one '@@' chiral specification in its isomeric SMILES (proxy for 5β configuration).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a 3-oxo-5β-steroid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned (helps in generating proper isomeric SMILES)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # --- Step 1: Identify fused, non-aromatic rings
    ring_info = mol.GetRingInfo()
    ring_atom_sets = []
    ring_sizes = []  # corresponding ring sizes
    # only use non-aromatic rings
    for ring in ring_info.AtomRings():
        if all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            ring_atom_sets.append(set(ring))
            ring_sizes.append(len(ring))
    if not ring_atom_sets:
        return False, "No non-aromatic rings found; no steroid nucleus"

    # Build a connectivity graph on rings: two rings are fused if they share at least 2 atoms.
    n = len(ring_atom_sets)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if len(ring_atom_sets[i].intersection(ring_atom_sets[j])) >= 2:
                graph[i].add(j)
                graph[j].add(i)
    
    # Identify connected fused ring systems (each a candidate for a steroid or related nucleus)
    visited = set()
    fused_components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    comp.add(current)
                    stack.extend(graph[current] - visited)
            fused_components.append(comp)
    
    candidate_nucleus_atoms = None
    nucleus_explanation = ""
    # For each fused component, try to find a subset of exactly 4 rings that satisfy the steroid nucleus signature.
    for comp in fused_components:
        if len(comp) < 4:
            continue  # not enough rings in this fused system
        comp_list = list(comp)
        # Consider all combinations of exactly 4 rings from this fused component.
        for subset in combinations(comp_list, 4):
            subset_sizes = [ring_sizes[i] for i in subset]
            if sorted(subset_sizes) != [5, 6, 6, 6]:
                continue
            # Form the union of atom indices for these four rings.
            atom_union = set()
            for i in subset:
                atom_union |= ring_atom_sets[i]
            carbon_count = sum(1 for idx in atom_union if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if carbon_count >= 17:
                candidate_nucleus_atoms = atom_union
                nucleus_explanation = (f"Fused nucleus found with rings of sizes {sorted(subset_sizes)} "
                                       f"and {carbon_count} carbons in the union of the candidate nucleus")
                break
        if candidate_nucleus_atoms is not None:
            break
    
    if candidate_nucleus_atoms is None:
        return False, "No fused ring system with a steroid nucleus signature (four rings with sizes [5,6,6,6] and ≥17 carbons) found"
    
    # --- Step 2: Check for a ring-bound ketone located within the candidate nucleus.
    # The SMARTS "[R]C(=O)[R]" matches a ketone with both neighbors being in rings.
    ketone_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_in_nucleus = False
    for match in ketone_matches:
        # match[0] is the carbonyl carbon.
        if match[0] in candidate_nucleus_atoms:
            ketone_in_nucleus = True
            break
    if not ketone_in_nucleus:
        return False, "No ring-bound ketone (3-oxo) group found within the candidate steroid nucleus"
    
    # --- Step 3: Check for beta configuration proxy using '@@' in the isomeric SMILES.
    iso_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    if "@@" not in iso_smi:
        return False, "No chiral center with '@@' (indicative of beta configuration) detected"
    
    explanation = ("Molecule has a candidate fused steroid nucleus (" + nucleus_explanation +
                   "), contains a ring-bound ketone (3-oxo) within that nucleus, and its isomeric SMILES '" +
                   iso_smi + "' contains '@@', indicative of 5β configuration.")
    return True, explanation

# Example usage (for testing, uncomment the following lines):
#if __name__ == "__main__":
#    tests = [
#        ("C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2C(=O)CO", "5beta-dihydrodeoxycorticosterone"),
#        # additional examples can be added here...
#    ]
#    for smi, name in tests:
#        result, reason = is_3_oxo_5beta_steroid(smi)
#        print(name, result, reason)