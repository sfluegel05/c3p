"""
Classifies: CHEBI:47787 11-oxo steroid
"""
#!/usr/bin/env python
"""
Classifies: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent at position 11.
This heuristic implementation first attempts to detect a steroid nucleus,
defined as four fused rings (three six-membered and one five-membered) having
a union of atoms in a slightly relaxed range. Next it verifies that one of the nucleus atoms
or one of its immediate neighbors bears a ketone (C=O) group.
Note:
  This method is heuristic and may miss exotic cases.
"""

from rdkit import Chem
from itertools import combinations

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 11-oxo steroid based on its SMILES string.
    It performs two steps:
      1) It searches for a steroid nucleus – defined as a set of 4 fused rings 
         (three six-membered and one five-membered) whose combined atoms number 
         is within a relaxed range (15 to 23 atoms).
      2) It checks whether a ketone group (C=O) is present either on the nucleus 
         or directly attached to an atom of the nucleus.
         
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 11-oxo steroid, False otherwise.
        str: A reason explaining the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get ring information and filter rings to sizes 5 or 6 (common in steroids).
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    candidate_rings = [set(r) for r in all_rings if len(r) in (5, 6)]
    
    if len(candidate_rings) < 4:
        return False, "Not enough 5- or 6-membered rings to form a steroid nucleus."
    
    # Look for a combination of 4 rings that have exactly one five-membered and three six-membered rings,
    # are mutually fused (each ring shares at least one atom with at least one other),
    # and the union of the atoms falls into a relaxed range (15 to 23 atoms).
    nucleus_atoms = None
    for ring_combo in combinations(candidate_rings, 4):
        count5 = sum(1 for r in ring_combo if len(r) == 5)
        count6 = sum(1 for r in ring_combo if len(r) == 6)
        if count5 != 1 or count6 != 3:
            continue  # Must have exactly one 5-membered and three 6-membered rings.
        # Check connectivity among these four rings:
        # Build a simple graph: nodes are rings; an edge exists if two rings share at least one atom.
        nodes = list(ring_combo)
        adj = {i: set() for i in range(4)}
        for i in range(4):
            for j in range(i+1, 4):
                if nodes[i] & nodes[j]:
                    adj[i].add(j)
                    adj[j].add(i)
        # Perform DFS to ensure all rings are connected.
        visited = set()
        stack = [0]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            stack.extend(adj[node] - visited)
        if len(visited) != 4:
            continue  # Not all rings are fused.
        # Compute the union of atoms from these rings.
        union_atoms = set()
        for r in ring_combo:
            union_atoms.update(r)
        if 15 <= len(union_atoms) <= 23:
            nucleus_atoms = union_atoms
            break  # Found a candidate steroid nucleus.
    
    if nucleus_atoms is None:
        return False, "No steroid nucleus (4 fused rings with appropriate sizes) found."
    
    # Prepare a SMARTS for a ketone group in a ring.
    # The pattern [#6;R][CX3](=O) looks for a ring-bound carbon attached to a carbonyl group.
    ketone_smarts = "[#6;R][CX3](=O)"
    ketone_pat = Chem.MolFromSmarts(ketone_smarts)
    if ketone_pat is None:
        return None, None  # SMARTS pattern creation failed.
    
    ketone_matches = mol.GetSubstructMatches(ketone_pat)
    if not ketone_matches:
        return False, "No ring-embedded ketone group detected (required for 11-oxo steroid)."
    
    # Check if any of the detected ketone groups are on the nucleus or directly attached to it.
    # Note: In the SMARTS our match, the second atom (index 1) is the carbonyl carbon.
    for match in ketone_matches:
        ketone_carbon = match[1]
        if ketone_carbon in nucleus_atoms:
            return True, "Steroid nucleus detected along with a ketone group on the nucleus (consistent with 11-oxo steroid)."
        # If the ketone carbon is not directly within the nucleus, check its neighbors.
        atom = mol.GetAtomWithIdx(ketone_carbon)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in nucleus_atoms:
                return True, ("Steroid nucleus detected and ketone group is immediately attached "
                              "to the nucleus (consistent with 11-oxo steroid).")
    
    return False, "Ketone group not located on or directly attached to the steroid nucleus (likely not 11-oxo steroid)."

# Example usage (you may uncomment to test):
# test_smiles = [
#     "O=C1C2=C([C@@]3(C(=O)C[C@@H]([C@]3(C1)C)C(=C)CCC(=O)O)C)[C@@H](O)C[C@@H]4[C@@]2(CCC(C4(C)C)=O)",  # 20(21)-dehydrolucidenic acid A
#     "[H][C@]12CC[C@@]3([H])[C@]4([H])CCC\\C(=C/CO)[C@@]4(C)CC(=O)[C@]3([H])[C@@]1(C)CC[C@@H](O)C2",  # 3alpha,21-dihydroxy-D-homo-5beta-pregn-17a(20)-en-11-one (previously false negative)
# ]
# for smi in test_smiles:
#     result, reason = is_11_oxo_steroid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")