"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D compounds (fat‐soluble hydroxy seco‐steroids).

Heuristic criteria (improved):
  1. SMILES must be valid.
  2. Molecular weight roughly within 250–700 Da.
  3. Must contain at least one hydroxyl group (-OH).
  4. Must contain a “seco‐steroid” core. We consider only rings of size 5 or 6 that
     are mostly carbocyclic (at least 4 of 5/6 atoms are carbons). Then, we build a graph
     of rings sharing atoms and require the largest connected (fused) component to have
     either 2 or 3 rings.
  5. Must contain at least one conjugated diene (C=C-C=C) pattern in which at least 2 atoms
     are not in any ring.
  6. Shows significant lipophilicity (calculated logP); we lower the threshold from 3.5 to 3.0
     to account for cases (e.g. glucuronides) that are slightly less lipophilic.
If a test fails, the function returns a corresponding reason.
Note: This is a heuristic filter and may misclassify edge cases.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound (fat‐soluble hydroxy seco‐steroid)
    based on its SMILES string.
    
    Heuristic criteria:
      - Valid SMILES.
      - Molecular weight roughly within 250 to 700 Da.
      - Contains at least one hydroxyl group (-OH).
      - Contains a secosteroid core defined as a fused ring system from rings of size 5 or 6
        that are largely carbocyclic (≥4 carbon atoms per ring) with a connected cluster of 2 or 3 rings.
      - Contains a conjugated diene (C=C-C=C) motif with at least two atoms not in any ring.
      - Has sufficient lipophilicity (calculated logP ≥ 3.0).
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule meets vitamin D heuristic criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Validate SMILES and create molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check molecular weight (≈ 250–700 Da).
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 700:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside the typical range (250–700 Da) for vitamin D compounds"
    
    # 3. Check for at least one hydroxyl group (-OH).
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "No hydroxyl (-OH) group found; vitamin D compounds must have at least one -OH"
    
    # 4. Identify candidate rings for the secosteroid core.
    # Consider only rings of size 5 or 6 that are “mostly” carbocyclic (≥4 carbons).
    ri = mol.GetRingInfo()
    all_rings = ri.AtomRings()  # each ring is a tuple of atom indices
    core_rings = []
    for ring in all_rings:
        if len(ring) in (5, 6):
            n_carbons = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
            if n_carbons >= 4:
                core_rings.append(set(ring))
    if not core_rings:
        return False, "No 5- or 6-membered (mostly carbocyclic) rings found; lacking a steroid-like core"
    
    # Build an adjacency graph: each ring is a node, and share an edge if they share >= 1 atom.
    n = len(core_rings)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if core_rings[i].intersection(core_rings[j]):
                adj[i].add(j)
                adj[j].add(i)
    
    # Find connected components (fused ring systems) using DFS.
    visited = set()
    components = []
    def dfs(node, comp):
        comp.add(node)
        for neigh in adj[node]:
            if neigh not in comp:
                dfs(neigh, comp)
    for i in range(n):
        if i not in visited:
            comp = set()
            dfs(i, comp)
            visited |= comp
            components.append(comp)
    
    # Find the largest component (by number of rings).
    largest_comp = max(components, key=lambda comp: len(comp))
    num_fused = len(largest_comp)
    if num_fused not in (2, 3):
        return False, f"Fused 5/6-membered ring system has {num_fused} rings; expected 2 or 3 for a vitamin D (seco-steroid) core"
    
    # 5. Look for a conjugated diene motif.
    # We will search using an atom-explicit pattern for a diene: four atoms with alternating double bonds.
    diene_smarts = Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]")
    diene_matches = mol.GetSubstructMatches(diene_smarts)
    found_diene = False
    for match in diene_matches:
        # Count how many of the atoms in the match are NOT in any ring.
        non_ring_atoms = sum(1 for idx in match if not mol.GetAtomWithIdx(idx).IsInRing())
        if non_ring_atoms >= 2:
            found_diene = True
            break
    if not found_diene:
        return False, "No conjugated diene (C=C-C=C) pattern found with sufficient exposure (≥ 2 atoms outside rings); may be missing the secosteroid signature"
    
    # 6. Calculate lipophilicity (logP) and check threshold.
    logP = Crippen.MolLogP(mol)
    if logP < 3.0:
        return False, f"Calculated logP ({logP:.1f}) is rather low; vitamin D compounds are expected to be lipophilic (logP ≥ 3.0)"
    
    return True, ("Molecule meets vitamin D criteria: appropriate molecular weight, contains at least one -OH, "
                   "a secosteroid core (2–3 fused 5/6-membered rings that are largely carbocyclic), a detectable conjugated diene motif, "
                   "and sufficient lipophilicity.")


# Example usage when running this module directly.
if __name__ == "__main__":
    # An example vitamin D compound (a vitamin D3 derivative)
    test_smiles = "C[C@H](C=C[C@H](C)C1CC[C@@H]2[C@@]1(CCCC2=CC=C3C[C@H](CCC3=C)O)C"
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)