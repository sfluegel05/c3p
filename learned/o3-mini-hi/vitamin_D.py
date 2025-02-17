"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D compounds (fat‐soluble hydroxy seco‐steroids).
This function uses simple heuristics:
  1. The input must be a valid SMILES string.
  2. The molecule’s exact molecular weight should lie in a typical range (250–600 Da).
  3. It must have at least one hydroxyl group (–OH).
  4. Many vitamin D compounds are “secosteroids” derived from cholesterol.
     Their core steroid nucleus (normally four fused rings) loses one ring 
     (usually the B ring); thus a fused ring analysis on rings of size 5 or 6 is used.
     If the largest connected component of such rings has exactly 3 fused rings,
     we accept the molecule as vitamin D.
If any of these conditions fail, the function gives a reason.
Note: Because structural definitions of vitamin D analogues are complex,
this algorithm gives a heuristic classification and may miss some edge cases.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Heuristic criteria:
      - Valid SMILES.
      - Molecular weight roughly within 250 to 600 Da.
      - Contains at least one hydroxyl group (-OH).
      - Contains a fused ring system (of 5– or 6-membered rings) forming 
        a secosteroid core with three fused rings, rather than the four rings 
        seen in typical full steroids.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule meets vitamin D heuristic criteria, False otherwise.
        str: Reason for classification.
    """
    # Convert SMILES to molecule and validate input.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (most vitamin D analogs cluster roughly in this range)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside the typical range for vitamin D"

    # Check that at least one hydroxyl group (-OH) is present.
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_smarts)
    if not hydroxyl_matches:
        return False, "No hydroxyl (-OH) group found; vitamin D compounds must have at least one -OH"

    # Get all ring systems from molecule
    # We'll only consider rings of size 5 or 6 (as seen in steroid cores).
    rings = mol.GetRingInfo().AtomRings()
    core_rings = [set(ring) for ring in rings if len(ring) in (5, 6)]
    if not core_rings:
        return False, "No 5- or 6-membered rings found; lacking steroid-like core"
    
    # Build an undirected graph between rings that share at least one atom.
    n = len(core_rings)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if core_rings[i].intersection(core_rings[j]):
                adj[i].add(j)
                adj[j].add(i)
    
    # Find the size of the largest connected component (fused ring system).
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
    
    # Vitamin D molecules are "seco-steroids" derived from full steroids.
    # Typical steroids have 4 fused rings; vitamin D compounds have lost one ring,
    # so we expect a fused system of about 3 rings in the steroid core.
    if largest != 3:
        return False, f"Fused 5/6-membered ring system has {largest} rings; expected 3 for a secosteroid (vitamin D) core"
    
    # (Optionally, one might also check fat-solubility; for instance, many vitamin D compounds have high logP.
    # Here we simply compute logP and check if it is greater than a threshold.)
    logP = Crippen.MolLogP(mol)
    if logP < 3:
        return False, f"Calculated logP ({logP:.1f}) is rather low; vitamin D compounds are expected to be lipophilic (high logP)"
    
    # If all heuristic tests pass, we consider it a vitamin D compound.
    return True, ("Molecule meets vitamin D criteria: appropriate molecular weight, has hydroxyl group(s), "
                  "exhibits a secosteroid core (3 fused 5/6-membered rings), and is lipophilic.")

# Example usage:
if __name__ == "__main__":
    test_smiles = "C[C@H](C=C[C@H](C)C1CC[C@@H]2[C@@]1(CCCC2=CC=C3C[C@H](CCC3=C)O)C"  # Example vitamin D3 SMILES
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)