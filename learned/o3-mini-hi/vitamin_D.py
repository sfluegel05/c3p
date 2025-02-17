"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D compounds (fat‐soluble hydroxy seco‐steroids).
This improved heuristic uses the following criteria:
  1. SMILES must be valid.
  2. Molecular weight roughly within 250–700 Da.
  3. Contains at least one hydroxyl group (-OH).
  4. Contains a “seco‐steroid” core defined here as the presence of a fused ring system
     consisting of 5- or 6-membered rings. We require that the largest connected component 
     (by sharing one or more atoms) contains 2 or 3 rings.
  5. Contains a conjugated diene (C=C-C=C) motif in which at least 2 of the 4 atoms do not lie 
     inside any ring, suggesting an open (seco) portion.
  6. Is lipophilic. We require a calculated logP greater than 3.5.
If any test fails, a reason is provided.
Note: This is a heuristic and might miss edge cases.
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
      - Contains a “seco‐steroid” core: a fused ring system (only counting rings of size 5 or 6)
        where the largest connected group has 2 or 3 rings.
      - Contains a conjugated diene (C=C-C=C) motif with at least two atoms not in any ring.
      - Has high lipophilicity (logP > 3.5).
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule meets the vitamin D heuristic criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Validate SMILES and create molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check molecular weight (250–700 Da).
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 700:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside the typical range (250–700) for vitamin D compounds"

    # 3. Check for at least one hydroxyl group (-OH).
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "No hydroxyl (-OH) group found; vitamin D compounds must have at least one -OH"

    # 4. Look for a steroid-like (seco-steroid) core.
    # Consider only rings of size 5 or 6.
    rings = mol.GetRingInfo().AtomRings()
    core_rings = [set(ring) for ring in rings if len(ring) in (5, 6)]
    if not core_rings:
        return False, "No 5- or 6-membered rings found; lacking a steroid-like core"

    # Build an adjacency (graph) of rings that share one or more atoms.
    n = len(core_rings)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if core_rings[i].intersection(core_rings[j]):
                adj[i].add(j)
                adj[j].add(i)

    # Use depth-first search (DFS) to find connected components among rings.
    visited = set()
    def dfs(node, component):
        component.add(node)
        for neigh in adj[node]:
            if neigh not in component:
                dfs(neigh, component)
                
    largest_component_size = 0
    for i in range(n):
        if i not in visited:
            component = set()
            dfs(i, component)
            visited.update(component)
            if len(component) > largest_component_size:
                largest_component_size = len(component)
    # Accept only if largest fused ring component is 2 or 3 rings.
    if largest_component_size not in (2, 3):
        return False, f"Fused 5/6-membered ring system has {largest_component_size} rings; expected 2 or 3 for a vitamin D (seco-steroid) core"

    # 5. Look for a conjugated diene pattern (C=C-C=C) outside ring system.
    diene_smarts = Chem.MolFromSmarts("C=C-C=C")
    diene_matches = mol.GetSubstructMatches(diene_smarts)
    found_diene = False
    for match in diene_matches:
        # Count the number of atoms in the match that are NOT in any ring.
        non_ring_atoms = sum(1 for idx in match if not mol.GetAtomWithIdx(idx).IsInRing())
        if non_ring_atoms >= 2:
            found_diene = True
            break
    if not found_diene:
        return False, "No conjugated diene (C=C-C=C) pattern found with sufficient exposure (at least 2 atoms out of rings); may be missing secosteroid signature"

    # 6. Check lipophilicity (logP).
    logP = Crippen.MolLogP(mol)
    if logP < 3.5:
        return False, f"Calculated logP ({logP:.1f}) is rather low; vitamin D compounds are expected to be lipophilic (logP > 3.5)"

    # All tests pass.
    return True, ("Molecule meets vitamin D criteria: appropriate molecular weight, at least one -OH, "
                   "a secosteroid core (2 or 3 fused 5/6-membered rings), the required conjugated diene motif, "
                   "and is lipophilic.")

# Example usage:
if __name__ == "__main__":
    # Testing with one vitamin D example.
    test_smiles = "C[C@H](C=C[C@H](C)C1CC[C@@H]2[C@@]1(CCCC2=CC=C3C[C@H](CCC3=C)O)C"
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)