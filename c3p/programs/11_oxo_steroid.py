"""
Classifies: CHEBI:47787 11-oxo steroid
"""
#!/usr/bin/env python
"""
Classifies: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent at position 11.
This heuristic implementation first attempts to detect a steroid nucleus,
defined as a fused system of 4 rings (three 6-membered and one 5-membered) with
a total union of about 15–20 atoms. Next it verifies that one of the nucleus atoms 
bears a ketone function (C=O) consistent with a 11-oxo substituent.
Note:
  This method is heuristic and may miss some steroids or mis‐classify exotic cases.
"""

from rdkit import Chem
from itertools import combinations

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 11-oxo steroid based on its SMILES string.
    It performs two steps:
      1) It searches for a steroid nucleus – defined as a set of 4 fused rings 
         (three six-membered and one five-membered) whose combined atoms number 
         is in an expected range (~15–20 atoms).
      2) It then checks whether there is a ketone group (C=O) embedded on one 
         of the atoms in that nucleus.

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
    
    # Obtain the ring information.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    # Filter rings: only consider rings of size 5 or 6 (common in steroids)
    candidate_rings = [set(r) for r in all_rings if len(r) in (5,6)]
    if len(candidate_rings) < 4:
        return False, "Not enough 5- or 6-membered rings to form a steroid nucleus."
    
    # We now try to find a combination of exactly 4 rings containing:
    #   - exactly one 5-membered ring and three 6-membered rings,
    #   - these rings must be fused (i.e. connected via shared atoms),
    #   - the union of atoms in these rings should be in the expected range (~15 to 20 atoms).
    nucleus_atoms = None
    for ring_combo in combinations(candidate_rings, 4):
        count5 = sum(1 for r in ring_combo if len(r)==5)
        count6 = sum(1 for r in ring_combo if len(r)==6)
        if count5 != 1 or count6 != 3:
            continue  # Must have exactly one 5-membered and three 6-membered rings.
        # Check connectivity among these four rings.
        # Build a simple graph: nodes = rings, and an edge exists if two rings share at least one atom.
        nodes = list(ring_combo)
        # Create an adjacency list for the 4 rings.
        adj = {i: set() for i in range(4)}
        for i in range(4):
            for j in range(i+1, 4):
                if nodes[i] & nodes[j]:
                    adj[i].add(j)
                    adj[j].add(i)
        # Use DFS to check connectivity.
        visited = set()
        stack = [0]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            stack.extend(adj[node] - visited)
        if len(visited) != 4:
            continue  # The four rings are not all fused.
        # Compute the union of atoms.
        union_atoms = set()
        for r in ring_combo:
            union_atoms.update(r)
        if 15 <= len(union_atoms) <= 21:
            nucleus_atoms = union_atoms
            break  # Found a valid steroid nucleus candidate.
    
    if nucleus_atoms is None:
        return False, "No steroid nucleus (4 fused rings with appropriate sizes) found."
    
    # Now, search for a ketone (C=O) group embedded in a ring.
    # We use a SMARTS that requires a carbonyl carbon (sp2) with a double-bonded oxygen,
    # and the carbon must be in a ring.
    ketone_smarts = "[#6;R][CX3](=O)"  # the first atom is ring-bound; the carbonyl (CX3) is our target.
    ketone_pat = Chem.MolFromSmarts(ketone_smarts)
    if ketone_pat is None:
        return None, None  # pattern creation failed.
    ketone_matches = mol.GetSubstructMatches(ketone_pat)
    if not ketone_matches:
        return False, "No ring-embedded ketone group detected (required for 11-oxo steroid)."
    
    # Check if any of the ketone groups involve an atom from the steroid nucleus.
    # In our SMARTS "[#6;R][CX3](=O)", the second atom (index 1) is the carbonyl carbon.
    for match in ketone_matches:
        ketone_carbon = match[1]
        if ketone_carbon in nucleus_atoms:
            return True, "Steroid nucleus detected along with a ketone group on the nucleus (consistent with 11-oxo steroid)."
    
    return False, "Ketone group not located on the steroid nucleus (likely not 11-oxo steroid)."

# Example usage (uncomment to test):
# smiles_list = [
#     "O=C1C2=C([C@@]3(C(=O)C[C@@H]([C@]3(C1)C)C(=C)CCC(=O)O)C)[C@@H](O)C[C@@H]4[C@@]2(CCC(C4(C)C)=O)C",  # 20(21)-dehydrolucidenic acid A
#     "C[C@]12CC(=O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@]2(O)C(=O)C=O",  # 21-dehydrocortisone
# ]
# for smi in smiles_list:
#     result, reason = is_11_oxo_steroid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")