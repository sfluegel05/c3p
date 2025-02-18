"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI: 3α-hydroxy steroid (improved heuristic)
A 3α-hydroxy steroid is defined as a steroid having a fused tetracyclic core (typically three six‐membered rings and one five‐membered ring)
with a hydroxyl (-OH) substituent on one of the six‐membered rings of that core.
The improved heuristic implemented here:
  1. Parses the SMILES and adds explicit hydrogens (to ensure –OH groups appear explicitly).
  2. Extracts all rings and keeps only those of size 5 or 6.
  3. Builds a connectivity graph among these rings (two rings are “fused” if they share at least two atoms).
  4. Finds the largest connected component (fusion group) and then requires that it consists of exactly 4 rings—
     with exactly three six-membered rings and one five-membered ring (the classical steroid nucleus).
  5. Finally, searches for any explicit –OH (using the [OX2H] SMARTS) that is attached to a carbon 
     (atomic number 6) from one of the six‐membered rings in that core.
Note: Determining a true “3α” orientation from SMILES stereochemistry is very challenging;
here we assume that having the –OH on the steroid core is a good enough proxy.
"""

from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3α-hydroxy steroid based on its SMILES string using an improved heuristic.
    
    The steps are:
      1. Parse the SMILES and add explicit hydrogens.
      2. Extract all rings and keep only those of size 5 or 6.
      3. Build a graph connecting rings that are fused (sharing >= 2 atoms) and identify the largest fused group.
      4. Verify that the fused ring system consists of exactly four rings: exactly three six-membered rings and one five-membered ring.
      5. Look for hydroxyl groups (SMARTS “[OX2H]”) and check that at least one –OH is attached to a carbon that is
         part of one of the six-membered rings in the fused steroid core.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3α-hydroxy steroid; False otherwise.
        str: A message explaining the classification decision.
    """
    # 1. Parse SMILES and add explicit hydrogens (to ensure –OH groups are visible)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # 2. Get ring information and select only rings of size 5 or 6.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    candidate_rings = [ring for ring in all_rings if len(ring) in (5, 6)]
    if not candidate_rings:
        return False, "No five or six-membered rings found (steroid core expected)"
    
    # 3. Build a graph where each candidate ring is a node.
    # Two rings are considered fused if they share at least 2 atoms.
    n = len(candidate_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        set_i = set(candidate_rings[i])
        for j in range(i+1, n):
            set_j = set(candidate_rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components (each is a fused ring group)
    visited = set()
    fused_groups = []  # will be list of sets (indices)
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
    
    if not fused_groups:
        return False, "No fused ring system found"
    
    # 4. Select the largest fused group and require that it consists of exactly four rings.
    largest_group = max(fused_groups, key=lambda comp: len(comp))
    if len(largest_group) != 4:
        return False, f"Fused ring system size ({len(largest_group)}) does not equal 4 rings (expected classical steroid core)"
    
    # Now, from the four candidate rings in the fused core, count how many are six-membered and five-membered.
    group_rings = [candidate_rings[i] for i in largest_group]
    six_count = 0
    five_count = 0
    six_ring_atoms = set()  # union of atom indices in six-membered rings
    for ring in group_rings:
        if len(ring) == 6:
            six_count += 1
            six_ring_atoms.update(ring)
        elif len(ring) == 5:
            five_count += 1
    # The classical steroid nucleus: three six-membered rings and one five-membered ring.
    if six_count != 3 or five_count != 1:
        return False, (f"Fused ring system does not match expected steroid core: "
                       f"found {six_count} six-membered rings and {five_count} five-membered rings")
    
    # 5. Identify hydroxyl groups.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    if oh_pattern is None:
        return False, "Error generating hydroxyl group pattern"
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if not oh_matches:
        return False, "No hydroxyl (-OH) groups found"
    
    # For each hydroxyl group, check if at least one is attached to a carbon atom that belongs to one of the six-membered rings
    # in the steroid core.
    for match in oh_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Check neighbors of the oxygen; look for a carbon (atomic num 6)
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue
            c_idx = neighbor.GetIdx()
            if c_idx in six_ring_atoms:
                # We found a hydroxyl attached to a carbon in one of the six-membered rings of the steroid core.
                return True, ("Molecule has the expected fused tetracyclic steroid core (3 six-membered and 1 five-membered rings) "
                              "and a hydroxyl group attached to a six-membered ring carbon, consistent with a 3α-hydroxy steroid")
    
    return False, "No hydroxyl group found on a six-membered ring carbon in the steroid core"

# Example usage:
if __name__ == "__main__":
    # Test on one of the provided examples: (3alpha,5alpha,17beta)-3-hydroxyandrostan-17-yl sulfate
    test_smiles = "C1[C@@]2([C@@]([C@@]3([C@](C[C@H](O)CC3)(C1)[H])C)(CC[C@@]4([C@H](CC[C@@]24[H])OS(O)(=O)=O)C)[H])[H]"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)