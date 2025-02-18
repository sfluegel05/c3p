"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
#!/usr/bin/env python3
"""
Classifies: 11β-hydroxy steroid
Defined as: Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta-configuration.
Heuristic approach:
  1. Parse the SMILES string.
  2. Get the molecule’s ring information and filter only rings of size 5 or 6,
     as these are the sizes found in a cyclopentanoperhydrophenanthrene nucleus.
  3. Build a fused ring graph (two rings are fused if they share at least 2 atoms)
     and group rings into connected components.
  4. For each component, require that:
       • it has at least 4 rings,
       • it contains at least 3 six-membered rings and at least one five-membered ring,
       • its union of ring atoms (when counting only carbons) is in the expected range (roughly 15–21 atoms).
  5. If a steroid‐like nucleus is found, then use the SMARTS "[C@@]([OX2H])" to detect beta–oriented –OH.
  6. Finally, require that at least one of these beta –OH groups is on a chiral carbon that is part
     of the steroid nucleus.
This should reduce false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is an 11β-hydroxy steroid based on its SMILES string.
    
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
    all_rings = list(mol.GetRingInfo().AtomRings())
    if not all_rings:
        return False, "No rings detected; not a steroid."
    
    # Filter rings: only consider rings of size 5 or 6 (steroid rings)
    valid_rings = [set(r) for r in all_rings if len(r) in (5,6)]
    if len(valid_rings) < 4:
        return False, f"Only {len(valid_rings)} rings of size 5 or 6 found; a steroid nucleus (>=4 rings) is expected."
    
    # Build fused graph: two rings are fused if they share at least 2 atoms.
    n = len(valid_rings)
    fused_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(valid_rings[i] & valid_rings[j]) >= 2:
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
    
    # Search for a steroid‐like nucleus in one of the components.
    steroid_cluster_atoms = None
    for comp in components:
        comp_rings = [valid_rings[i] for i in comp]
        if len(comp_rings) < 4:
            continue
        
        # Count ring sizes.
        count6 = sum(1 for ring in comp_rings if len(ring) == 6)
        count5 = sum(1 for ring in comp_rings if len(ring) == 5)
        
        if count6 < 3 or count5 < 1:
            continue
        
        # Combine all atoms (the nucleus)
        comp_atoms = set()
        for ring in comp_rings:
            comp_atoms |= ring
        # Check that the nucleus has a reasonable size:
        # Steroid core usually has about 17 carbon atoms. Allow some wiggle room.
        carbon_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
        nucleus_carbons = comp_atoms & set(carbon_atoms)
        if not (15 <= len(nucleus_carbons) <= 21):
            continue
        
        steroid_cluster_atoms = comp_atoms  # record the atom indices involved in the nucleus
        break
    
    if steroid_cluster_atoms is None:
        return False, "No fused steroid nucleus (with >=3 six-membered rings and 1 five-membered ring and appropriate carbon count) detected."
    
    # Look for beta-oriented -OH using SMARTS.
    # [C@@] means chiral carbon with beta configuration; ([OX2H]) means directly bonded -OH.
    beta_oh_smarts = "[C@@]([OX2H])"
    beta_oh_query = Chem.MolFromSmarts(beta_oh_smarts)
    if beta_oh_query is None:
        return False, "Error constructing SMARTS query for beta-oriented hydroxyl group."
    
    matches = mol.GetSubstructMatches(beta_oh_query)
    if not matches:
        return False, "No beta-oriented hydroxyl group (SMARTS '[C@@]([OX2H])') found in the molecule."
    
    # Check that at least one match is on an atom that is in the steroid nucleus.
    for match in matches:
        # match[0] is the chiral carbon
        if match[0] in steroid_cluster_atoms:
            return True, "Steroid nucleus detected (fused rings with appropriate sizes and carbon count) and beta-oriented hydroxyl group found on it."
    
    return False, "Beta-oriented hydroxyl group not located on the identified steroid nucleus."

# Example usage and tests.
if __name__ == "__main__":
    # A few test SMILES examples (including some from the provided list)
    test_smiles_list = [
        # prednisolone phosphate
        "C[C@@]12C[C@@H]([C@]3([C@]([C@]2([H])CC[C@@]1(C(COP(=O)(O)O)=O)O)(CCC=4[C@]3(C)C=CC(C4)=O)[H])[H])O",
        # 5beta-dihydroaldosterone
        "C[C@]12CCC(=O)C[C@H]1CC[C@H]1[C@@H]3CC[C@H](C(=O)CO)[C@]3(C[C@H](O)[C@H]21)C=O",
        # withalongolide C
        "C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3[C@@H](O)C[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1",
        # A false positive example (a carbohydrate or non-steroid fused ring system)
        "CO[C@H]1C[C@H](O[C@H]2CC[C@@]3(C)[C@H](CC[C@@H]4[C@@H]3[C@H](O)C[C@]3(C)[C@H](CC[C@]43O)C3=CC(=O)OC3)C2)O[C@@H](C)[C@@H]1O"
    ]
    
    for sm in test_smiles_list:
        result, reason = is_11beta_hydroxy_steroid(sm)
        print("SMILES:", sm)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)