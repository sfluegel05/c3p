"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI: 3α-hydroxy steroid (improved heuristic)
A 3α-hydroxy steroid is defined as a steroid having a fused tetracyclic core (typically three six‐membered rings and one five‐membered ring)
with a hydroxyl (-OH) substituent on one of the six‐membered rings of that core.
This improved heuristic:
  1. Adds explicit hydrogens so that hydroxyl groups are correctly recognized.
  2. Identifies only rings of size 5 or 6.
  3. Builds a connectivity graph among these rings (an edge exists when two rings share at least 2 atoms) to identify a fused core.
  4. Checks that the largest fused group contains at least 4 rings (and, typically, three 6‑membered rings and one 5‑membered ring).
  5. Then ensures that at least one hydroxyl group (–OH) is attached to a carbon that is a member of a six‐membered ring in the core.
"""

from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3α-hydroxy steroid based on its SMILES string using an improved heuristic.
    
    The steps are:
      1. Parse the SMILES and add explicit hydrogens.
      2. Extract ring information, but only keep rings of size 5 or 6.
      3. Build a graph connecting rings that are fused (sharing at least 2 atoms) and select the largest connected component.
      4. Verify the fused core is consistent with a steroid (>=4 rings; typically >=3 rings of size 6 and >=1 ring of size 5).
      5. Find hydroxyl groups via the “[OX2H]” SMARTS and check if any is attached to a carbon that belongs to a six-membered ring
         of the fused core.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3α-hydroxy steroid; False otherwise.
        str: A message explaining the reasoning behind the classification.
    """
    # Parse and add hydrogens to ensure –OH groups are explicit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    # Keep only rings of size 5 or 6 as candidate steroid rings.
    candidate_rings = [ring for ring in all_rings if len(ring) in (5,6)]
    if not candidate_rings:
        return False, "No five or six-membered rings found (steroid core expected)"
    
    # Build a graph among candidate rings.
    # Two rings are considered fused if they share at least two atoms.
    n = len(candidate_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        set_i = set(candidate_rings[i])
        for j in range(i+1, n):
            set_j = set(candidate_rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components (fused groups) in the ring graph.
    visited = set()
    fused_groups = []
    for i in range(n):
        if i in visited:
            continue
        comp = set()
        stack = [i]
        while stack:
            cur = stack.pop()
            if cur in comp:
                continue
            comp.add(cur)
            for neighbor in ring_graph[cur]:
                if neighbor not in comp:
                    stack.append(neighbor)
        visited |= comp
        fused_groups.append(comp)
    
    # Select the largest fused group
    if not fused_groups:
        return False, "No fused ring system found"
    largest_group = max(fused_groups, key=lambda comp: len(comp))
    
    # Count rings in the largest fused group and remember six-membered ring atom indices.
    group_rings = [candidate_rings[i] for i in largest_group]
    num_rings = len(group_rings)
    six_ring_atoms = set()
    six_count = 0
    five_count = 0
    for ring in group_rings:
        if len(ring) == 6:
            six_ring_atoms.update(ring)
            six_count += 1
        elif len(ring) == 5:
            five_count += 1

    # A steroid should have a fused system with at least 4 rings.
    if num_rings < 4:
        return False, f"Fused ring system too small: {num_rings} rings (expected at least 4 in a steroid core)"
    # Typical steroid core has 3 six-membered rings and 1 five-membered ring.
    if six_count < 3 or five_count < 1:
        return False, (f"Fused ring system does not match the expected steroid core pattern "
                       f"(found {six_count} six-membered and {five_count} five-membered rings)")
    
    # Now search for hydroxyl groups.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    if oh_pattern is None:
        return False, "Error generating hydroxyl pattern"
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if not oh_matches:
        return False, "No hydroxyl (-OH) groups found in the molecule"
    
    # For each hydroxyl group, check if its attached carbon is part of a six-membered ring in the fused core.
    for match in oh_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue
            c_idx = neighbor.GetIdx()
            if c_idx in six_ring_atoms:
                return True, ("Molecule has a fused tetracyclic steroid core ("
                              f"{num_rings} fused rings: {six_count} six-membered and {five_count} five-membered) "
                              "and contains a hydroxyl group attached to a six-membered ring carbon in the core, "
                              "consistent with a 3α-hydroxy steroid")
    
    return False, "No hydroxyl group is attached to a six-membered ring carbon in the fused steroid core"


# Example usage:
if __name__ == "__main__":
    # Test on one of the provided examples: (3alpha,5alpha,17beta)-3-hydroxyandrostan-17-yl sulfate
    test_smiles = "C1[C@@]2([C@@]([C@@]3([C@](C[C@H](O)CC3)(C1)[H])C)(CC[C@@]4([C@H](CC[C@@]24[H])OS(O)(=O)=O)C)[H])[H]"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)