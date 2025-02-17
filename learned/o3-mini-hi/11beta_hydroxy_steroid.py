"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
#!/usr/bin/env python3
"""
Classifies: 11β-hydroxy steroid
Defined as: Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta- configuration.

Heuristic improvements:
  1. Identify fused ring clusters. We build a graph where each ring (taken from mol.GetRingInfo())
     is a node and rings are connected (“fused”) if they share at least 2 atoms.
  2. Among each connected component, we require that it contains at least 4 rings with at least
     3 six-membered rings and at least 1 five-membered ring (typical of a steroid nucleus,
     i.e. the cyclopentanoperhydrophenanthrene skeleton).
  3. Then we look for a beta-oriented hydroxyl group. We use the SMARTS "[C@@]([OX2H])" to select
     chiral carbons with the beta (“@@”) configuration bonded directly to a –OH group.
  4. Finally, we insist that at least one match from (3) exists on an atom that belongs to the fused
     steroid nucleus.
     
Note: This is a heuristic approach that may not be perfect.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is an 11β-hydroxy steroid based on its SMILES string.
    
    Heuristics:
      1. Identify fused ring clusters from the molecule's ring information. For a steroid nucleus,
         we require at least one connected set of rings that has:
           - At least 4 rings total,
           - At least 3 rings of size 6 and 
           - At least 1 ring of size 5 (as in the cyclopentanoperhydrophenanthrene group).
      2. Look for a beta-oriented hydroxyl group using the SMARTS "[C@@]([OX2H])" pattern.
      3. Check that at least one such beta –OH chiral carbon is part of the fused steroid nucleus.
      
    Returns:
      (True, reason) if the criteria are met,
      (False, reason) if not.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string."
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"
    
    # Get all ring information (each ring is a tuple of atom indices)
    rings = list(mol.GetRingInfo().AtomRings())
    if not rings or len(rings) < 4:
        return False, f"Found only {len(rings)} rings; steroid nucleus (fused 4 rings) expected."
    
    # Create a list of sets for easier membership checking
    ring_sets = [set(r) for r in rings]
    
    # Build a graph of rings: two rings are "fused" if they share at least 2 atoms.
    n = len(ring_sets)
    fused_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(ring_sets[i] & ring_sets[j]) >= 2:
                fused_graph[i].add(j)
                fused_graph[j].add(i)
                
    # Find connected components among the rings.
    visited = set()
    components = []
    
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    comp.add(current)
                    stack.extend(fused_graph[current] - visited)
            components.append(comp)

    # Look for a steroid-like fused ring cluster.
    steroid_cluster_atoms = set()
    cluster_found = False
    for comp in components:
        comp_rings = [ring_sets[i] for i in comp]
        num_rings = len(comp_rings)
        count6 = sum(1 for ring in comp_rings if len(ring) == 6)
        count5 = sum(1 for ring in comp_rings if len(ring) == 5)
        # Typical steroid nucleus: 4 fused rings (with 3 six-membered and 1 five-membered ring)
        if num_rings >= 4 and count6 >= 3 and count5 >= 1:
            # Combine all atoms in these rings
            for r in comp_rings:
                steroid_cluster_atoms |= r
            cluster_found = True
            break
            
    if not cluster_found:
        return False, "No fused steroid nucleus (with ≥3 six-membered rings & 1 five-membered ring) detected."
    
    # Look for beta-oriented -OH using SMARTS.
    # [C@@] means a chiral carbon with the beta descriptor and ([OX2H]) means directly bonded to an -OH.
    beta_oh_smarts = "[C@@]([OX2H])"
    beta_oh_query = Chem.MolFromSmarts(beta_oh_smarts)
    if beta_oh_query is None:
        return False, "Error constructing SMARTS query for beta-oriented hydroxyl group."
    
    matches = mol.GetSubstructMatches(beta_oh_query)
    if not matches:
        return False, "No beta-oriented hydroxyl group (SMARTS '[C@@]([OX2H])') found in the molecule."
    
    # Ensure that at least one match is on an atom within the fused steroid nucleus.
    for match in matches:
        carbon_idx = match[0]  # first atom (the chiral carbon) in the match
        if carbon_idx in steroid_cluster_atoms:
            return True, "Steroid nucleus detected (fused cluster with requisite rings) and beta-oriented hydroxyl group found on it."
    
    return False, "Beta-oriented hydroxyl group not found on the steroid nucleus."

# Example usage and tests.
if __name__ == "__main__":
    # A few test SMILES examples (you can add more)
    test_smiles_list = [
        # prednisolone phosphate
        "C[C@@]12C[C@@H]([C@]3([C@]([C@]2([H])CC[C@@]1(C(COP(=O)(O)O)=O)O)(CCC=4[C@]3(C)C=CC(C4)=O)[H])[H])O",
        # 5beta-dihydroaldosterone
        "C[C@]12CCC(=O)C[C@H]1CC[C@H]1[C@@H]3CC[C@H](C(=O)CO)[C@]3(C[C@H](O)[C@H]21)C=O",
        # withalongolide C
        "C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3[C@@H](O)C[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1",
        # a false positive example (a carbohydrate, not a steroid)
        "O([C@H]1[C@H](O)[C@H](O[C@H](O)[C@H]1O)CO)[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)CO"
    ]
    
    for sm in test_smiles_list:
        result, reason = is_11beta_hydroxy_steroid(sm)
        print("SMILES:", sm)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)