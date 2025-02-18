"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3β-hydroxy steroid
Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β‐position.
The program identifies a steroid nucleus by looking for a fused tetracyclic ring system
composed of one 5-membered and three 6-membered rings (each predominantly carbon)
that are connected (fused by sharing at least two atoms). Then, it checks that at least one
beta‐oriented hydroxyl group ([C@@H](O)) is attached to an atom in that system.
"""

from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.
    The classifier verifies:
      1. The presence of a steroid nucleus (i.e. a fused tetracyclic ring system with one 5-ring and three 6-rings,
         where the rings are largely carbon).
      2. The presence of at least one beta-oriented hydroxyl group ([C@@H](O)) attached to an atom within the steroid nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3β-hydroxy steroid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information as a list of tuples (each tuple holds atom indices for a ring)
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Filter rings based on size and composition.
    # We expect rings in steroids to be 5- or 6-membered and mostly carbons.
    filtered_rings = []
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        # Count the carbons in the ring
        n_carbon = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        # If at least 80% of the atoms are carbons, keep this ring.
        if n_carbon / len(ring) >= 0.8:
            filtered_rings.append(ring)
    
    if not filtered_rings:
        return False, "No appropriate rings found for steroid nucleus"
    
    # Build a graph of rings where each ring is a node.
    # Two rings are considered fused if they share at least two atoms.
    ring_graph = {i: set() for i in range(len(filtered_rings))}
    for i in range(len(filtered_rings)):
        set_i = set(filtered_rings[i])
        for j in range(i+1, len(filtered_rings)):
            set_j = set(filtered_rings[j])
            if len(set_i & set_j) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Now find connected components in the ring graph.
    seen = set()
    fused_components = []
    for i in ring_graph:
        if i in seen:
            continue
        stack = [i]
        component = set()
        while stack:
            cur = stack.pop()
            if cur in component:
                continue
            component.add(cur)
            for neighbor in ring_graph[cur]:
                if neighbor not in component:
                    stack.append(neighbor)
        seen |= component
        fused_components.append(component)
    
    # Look for a fused component that has at least 4 rings.
    steroid_component = None
    for comp in fused_components:
        if len(comp) >= 4:
            # Count how many are 5-membered and how many are 6-membered
            n5 = 0
            n6 = 0
            for idx in comp:
                if len(filtered_rings[idx]) == 5:
                    n5 += 1
                elif len(filtered_rings[idx]) == 6:
                    n6 += 1
            if n5 >= 1 and n6 >= 3:
                steroid_component = comp
                break
    if steroid_component is None:
        return False, "Steroid nucleus not found (no fused tetracyclic system with expected ring sizes)"
    
    # Combine all atom indices belonging to the steroid nucleus.
    steroid_atoms = set()
    for idx in steroid_component:
        steroid_atoms.update(filtered_rings[idx])
    
    # Define SMARTS for a beta-oriented hydroxyl group.
    # This typically is encoded as [C@@H](O)
    beta_hydroxy_smarts = "[C@@H](O)"
    beta_hydroxy_query = Chem.MolFromSmarts(beta_hydroxy_smarts)
    
    # Find matches for beta-hydroxy group in the molecule.
    beta_matches = mol.GetSubstructMatches(beta_hydroxy_query)
    if not beta_matches:
        return False, "Beta-oriented hydroxyl group not found"
    
    # Require that at least one of these beta-hydroxy carbons is part of the steroid nucleus.
    found_beta_in_nucleus = False
    for match in beta_matches:
        # match is a tuple of atom indices corresponding to the query pattern;
        # the first atom should be the carbon with the OH.
        if match[0] in steroid_atoms:
            found_beta_in_nucleus = True
            break
    if not found_beta_in_nucleus:
        return False, "Beta-oriented hydroxyl group not attached to the steroid nucleus"
    
    return True, "Molecule contains a steroid nucleus with a beta-oriented (3β) hydroxyl group"

# Example usage:
# test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@@H](O)C(C)C"  # (24R)-24-hydroxycholesterol
# result, reason = is_3beta_hydroxy_steroid(test_smiles)
# print(result, reason)