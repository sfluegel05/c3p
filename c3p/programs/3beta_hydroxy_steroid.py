"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3β-hydroxy steroid
Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β‐position.

Our improved classifier proceeds as follows:
  1. Parse the molecule and obtain rings.
  2. Filter to rings of sizes 5 or 6 in which at least 67% of the atoms are carbon.
  3. Build a fused ring graph (where rings sharing at least 2 atoms are connected).
  4. For each connected fused ring component, require:
       - At least 4 rings are present.
       - At least one of the rings is 5-membered and at least three are 6-membered.
       - The union of atoms (the candidate steroid nucleus) is of a reasonable size (15–23 atoms)
         and is predominantly (>=80%) carbon.
  5. Look for at least one beta–oriented hydroxyl group ([C@@H](O)) whose carbon is within the nucleus.
  
If all conditions are met, we classify the molecule as a 3β–hydroxy steroid.
"""

from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.

    The classifier checks:
      1. That a steroid nucleus is present – by finding a fused tetracyclic system (4 or more rings)
         with at least one 5-membered and three 6-membered rings. In addition, the union of atoms 
         making up the cluster (“nucleus”) must be of a size consistent with steroids (15–23 atoms)
         and mostly carbon (>=80%).
      2. That at least one beta–oriented hydroxyl group ([C@@H](O)) is attached to an atom of the nucleus.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3β–hydroxy steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain all rings from the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Step 1: Filter for rings of size 5 or 6 with at least 67% carbon.
    filtered_rings = []
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbons / len(ring) >= 0.67:
            filtered_rings.append(ring)
    if not filtered_rings:
        return False, "No rings of size 5 or 6 with sufficient carbon content found"
    
    # Step 2: Build graph of fused rings (rings sharing at least 2 atoms are fused).
    ring_graph = {i: set() for i in range(len(filtered_rings))}
    for i in range(len(filtered_rings)):
        ring_i = set(filtered_rings[i])
        for j in range(i+1, len(filtered_rings)):
            ring_j = set(filtered_rings[j])
            if len(ring_i & ring_j) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
                
    # Gather connected components from the fused ring graph.
    seen = set()
    fused_components = []
    for i in ring_graph:
        if i in seen:
            continue
        stack = [i]
        comp = set()
        while stack:
            node = stack.pop()
            if node in comp:
                continue
            comp.add(node)
            for neighbor in ring_graph[node]:
                if neighbor not in comp:
                    stack.append(neighbor)
        seen |= comp
        fused_components.append(comp)
    
    # Step 3: Look for a fused component with at least 4 rings,
    # with at least one 5-membered ring and at least three 6-membered rings.
    # Also require that the nucleus (the union of ring atoms) is of a reasonable size (15-23 atoms)
    # and that at least 80% of its atoms are carbons.
    steroid_nucleus = None
    for comp in fused_components:
        if len(comp) < 4:
            continue
        n5 = sum(1 for idx in comp if len(filtered_rings[idx]) == 5)
        n6 = sum(1 for idx in comp if len(filtered_rings[idx]) == 6)
        if n5 < 1 or n6 < 3:
            continue

        # Get the union of atoms corresponding to this fused component.
        nucleus_atoms = set()
        for idx in comp:
            nucleus_atoms.update(filtered_rings[idx])
            
        # Check the nucleus size is in a typical steroid core range.
        if not (15 <= len(nucleus_atoms) <= 23):
            # If the fused component is too big, it might be part of a decorated ring system.
            continue
            
        # Check that the nucleus is predominantly carbon.
        nucleus_carbons = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if nucleus_carbons / len(nucleus_atoms) < 0.80:
            continue
            
        steroid_nucleus = nucleus_atoms
        break

    if steroid_nucleus is None:
        return False, "Steroid nucleus not found (no appropriate fused tetracyclic cluster detected)"

    # (Optional) Also check that the nucleus constitutes a reasonable fraction of the molecule.
    total_atoms = mol.GetNumAtoms()
    if len(steroid_nucleus)/total_atoms < 0.20:
        return False, "Steroid nucleus constitutes less than 20% of the molecule's atoms"

    # Step 4: Search for beta-oriented hydroxyl groups.
    # Use SMARTS "[C@@H](O)" which defines a chiral carbon with attached hydroxyl.
    beta_oh_query = Chem.MolFromSmarts("[C@@H](O)")
    beta_matches = mol.GetSubstructMatches(beta_oh_query)
    if not beta_matches:
        return False, "No beta-oriented hydroxyl group ([C@@H](O)) found"
    
    # Verify that at least one beta-OH is attached to an atom of the steroid nucleus.
    for match in beta_matches:
        # match[0] is the carbon with defined stereochemistry.
        if match[0] in steroid_nucleus:
            return True, "Molecule contains a steroid nucleus with a beta-oriented (3β) hydroxyl group"
    
    return False, "Beta–oriented hydroxyl group not attached to the steroid nucleus"

# Example usage:
# test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@@H](O)C(C)C"
# result, reason = is_3beta_hydroxy_steroid(test_smiles)
# print(result, reason)