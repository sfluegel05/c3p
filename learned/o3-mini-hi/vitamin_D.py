"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D compounds (fat‐soluble hydroxy seco‐steroids).
This improved heuristic uses the following criteria:
  1. SMILES must be valid.
  2. Molecular weight within 250–700 Da.
  3. At least one hydroxyl group (-OH) is present.
  4. The molecule shows evidence of a “seco‐steroid” core – a partial fused ring system 
     of 5– or 6-membered rings. We allow the largest connected component (by shared atoms)
     to have either 2 or 3 rings.
  5. The compound must have at least one conjugated diene pattern (C=C–C=C) that is not completely 
     inside a ring, suggesting the open (seco) ring.
  6. The molecule is lipophilic, with a calculated logP greater than 4.
If any test fails, a reason is provided.
Note: This is a heuristic approach and may miss edge cases.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Heuristic criteria:
      - Valid SMILES.
      - Molecular weight roughly within 250 to 700 Da.
      - Contains at least one hydroxyl group (-OH).
      - Contains a secosteroid (vitamin D) core: a fused ring system (of 5- or 6-membered rings)
        with 2 or 3 connected rings.
      - Contains a conjugated diene pattern (C=C-C=C) with at least one atom not in a ring.
      - Has high lipophilicity (logP > 4).
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule meets the vitamin D heuristic criteria, False otherwise.
        str: Reason for classification or the failure reason.
    """
    # 1. Convert SMILES to molecule and validate input.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 700:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside the typical range (250–700) for vitamin D compounds"

    # 3. Check that at least one hydroxyl group (-OH) is present.
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "No hydroxyl (-OH) group found; vitamin D compounds must have at least one -OH"

    # 4. Analyze fused ring systems (only consider rings of size 5 or 6).
    rings = mol.GetRingInfo().AtomRings()
    core_rings = [set(ring) for ring in rings if len(ring) in (5, 6)]
    if not core_rings:
        return False, "No 5- or 6-membered rings found; lacking a steroid-like core"
    
    # Build a graph of rings that share at least one atom.
    n = len(core_rings)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if core_rings[i].intersection(core_rings[j]):
                adj[i].add(j)
                adj[j].add(i)
    
    # DFS to find largest connected component of rings.
    visited = set()
    def dfs(node, comp):
        comp.add(node)
        for neigh in adj[node]:
            if neigh not in comp:
                dfs(neigh, comp)
    
    largest = 0
    for i in range(n):
        if i not in visited:
            comp = set()
            dfs(i, comp)
            visited.update(comp)
            if len(comp) > largest:
                largest = len(comp)
    # Accept if the fused ring count is 2 or 3.
    if largest not in (2, 3):
        return False, f"Fused 5/6-membered ring system has {largest} rings; expected 2 or 3 for a vitamin D (seco-steroid) core"

    # 5. Look for a conjugated diene (C=C-C=C) pattern.
    # This pattern is indicative of the broken ring found in secosteroids.
    conjugated_diene_smarts = Chem.MolFromSmarts("C=C-C=C")
    diene_matches = mol.GetSubstructMatches(conjugated_diene_smarts)
    found_conjugated_diene = False
    for match in diene_matches:
        # Check if at least one atom in the match is not in any ring.
        if any(not mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            found_conjugated_diene = True
            break
    if not found_conjugated_diene:
        return False, "No conjugated diene (C=C-C=C) pattern found outside rings; might be missing the secosteroid signature"

    # 6. Check lipophilicity.
    logP = Crippen.MolLogP(mol)
    if logP < 4:
        return False, f"Calculated logP ({logP:.1f}) is rather low; vitamin D compounds are expected to be lipophilic (logP > 4)"

    # If all tests pass, accept the compound.
    return True, ("Molecule meets vitamin D criteria: appropriate molecular weight, at least one -OH, "
                  "a secosteroid core (2 or 3 fused 5/6-membered rings), a conjugated diene motif, and is lipophilic.")

# Example usage:
if __name__ == "__main__":
    # Example vitamin D compound SMILES (one of those provided):
    test_smiles = "C[C@H](C=C[C@H](C)C1CC[C@@H]2[C@@]1(CCCC2=CC=C3C[C@H](CCC3=C)O)C"
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)