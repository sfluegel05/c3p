"""
Classifies: CHEBI:47788 3-oxo steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo steroid 
Definition: Any oxo steroid where an oxo (ketone) substituent is located at position 3.
This script uses a heuristic:
  1. The molecule must have a fused ring system typical of steroids (a set of 4 rings – usually three 6-membered and one 5-membered – that share bonds).
  2. At least one ketone group ([#6][CX3](=O)[#6]) must be present and the ketone C atom must belong to one of the six-membered rings of the fused steroid system.
Note: This is an approximate and heuristic method.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is defined (heuristically) as a molecule with a steroid-like fused
    tetracyclic ring system AND having a ketone group (C(=O)) located on a 6-membered ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a 3-oxo steroid, False otherwise.
        str : Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Obtain ring information (each ring is given as a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found, not a steroid"

    # Filter rings to only include those of size 5 or 6
    candidate_rings = [ring for ring in ring_info if len(ring) in (5,6)]
    if len(candidate_rings) < 4:
        return False, "Less than 4 rings of size 5 or 6, not a typical steroid nucleus"

    # Build a 'ring graph': nodes are ring indices and an edge between rings that share 2 or more atoms (i.e. fused)
    n = len(candidate_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # if rings i and j share at least 2 atoms, they are fused
            if len(set(candidate_rings[i]).intersection(candidate_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find the largest connected component of fused rings.
    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                current = stack.pop()
                if current not in comp:
                    comp.add(current)
                    visited.add(current)
                    stack.extend(ring_graph[current] - comp)
            components.append(comp)
    # For a steroid, we typically expect a fused system of 4 rings.
    largest_comp = max(components, key=lambda c: len(c))
    if len(largest_comp) < 4:
        return False, "No fused ring system of 4 rings detected, not a steroid nucleus"
    
    # Collect all atom indices that belong to six-membered rings in this largest fused set.
    six_membered_atoms = set()
    for idx in largest_comp:
        ring = candidate_rings[idx]
        if len(ring) == 6:
            six_membered_atoms.update(ring)

    if not six_membered_atoms:
        return False, "No six-membered rings found in the fused system, cannot locate typical 3-oxo site"

    # Use SMARTS to find ketone groups: a carbon bonded to O via a double bond
    # The pattern "[#6][CX3](=O)[#6]" ensures that the carbonyl carbon (the [CX3](=O) part)
    # is flanked by carbons.
    ketone_smarts = "[#6][CX3](=O)[#6]"
    ketone_query = Chem.MolFromSmarts(ketone_smarts)
    ketone_matches = mol.GetSubstructMatches(ketone_query)
    if not ketone_matches:
        return False, "No ketone groups (C(=O) flanked by carbons) found"

    # For each ketone match, the middle atom (index 1 in the match tuple) is the carbonyl carbon.
    # If one of those carbon atoms is in a six-membered ring of the fused system, we infer a likely 3-oxo steroid.
    for match in ketone_matches:
        carbonyl_idx = match[1]
        if carbonyl_idx in six_membered_atoms:
            return True, "Found ketone group on a six-membered ring of the steroid nucleus (likely at position 3)"
    
    return False, "Ketone group not located on a six-membered ring of the steroid nucleus, not a 3-oxo steroid"

# Example usage (uncomment for testing):
# smiles_list = [
#    "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO",  # methyl prednisolone-16alpha-carboxylate
#    "C[C@]12CC[C@H]3[C@@H](CC=C4CC(=O)CC[C@]34C)[C@@H]1CCC2=O",  # androst-5-ene-3,17-dione
# ]
# for smi in smiles_list:
#     result, reason = is_3_oxo_steroid(smi)
#     print(smi[:50], "->", result, "|", reason)