"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3β-hydroxy steroid
Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β‐position.
The program first looks for a steroid nucleus – defined as a fused tetracyclic ring system 
composed of one 5-membered and three 6-membered rings (with the rings being predominantly carbon) 
– and then it confirms that at least one beta‐oriented hydroxyl group ([C@@H](O)) is attached 
to an atom within that nucleus. Additionally, to avoid picking up large conjugated steroids 
(which are not our target simple 3β-hydroxy steroids), the code verifies that the nucleus constitutes 
at least 50% of the atoms in the molecule.
"""

from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.
    The classifier checks:
      1. Whether a steroid nucleus is present (a fused tetracyclic system with one 5-membered and three 6-membered rings).
         (Rings are filtered to have either 5 or 6 atoms and at least 67% carbons.)
      2. That the nucleus is not “diluted” by many extra substituents (the atoms in the nucleus must represent at least 50% of the molecule).
      3. That there is at least one beta-oriented hydroxyl group ([C@@H](O)) attached to an atom of the nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3β-hydroxy steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings (as tuples of atom indices) from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Filter rings: allow only rings of size 5 or 6 with at least ~67% carbons.
    filtered_rings = []
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        n_carbon = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbon / len(ring) >= 0.67:
            filtered_rings.append(ring)
    
    if not filtered_rings:
        return False, "No appropriate rings found for steroid nucleus"
    
    # Build a graph of rings. Two rings are “fused” if they share at least 2 atoms.
    ring_graph = {i: set() for i in range(len(filtered_rings))}
    for i in range(len(filtered_rings)):
        set_i = set(filtered_rings[i])
        for j in range(i+1, len(filtered_rings)):
            set_j = set(filtered_rings[j])
            if len(set_i & set_j) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components (fused ring clusters) in the ring graph.
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
            # Count how many 5-membered and 6-membered rings are in this component.
            n5 = sum(1 for idx in comp if len(filtered_rings[idx]) == 5)
            n6 = sum(1 for idx in comp if len(filtered_rings[idx]) == 6)
            if n5 >= 1 and n6 >= 3:
                steroid_component = comp
                break
    if steroid_component is None:
        return False, "Steroid nucleus not found (no fused tetracyclic system with expected ring sizes)"
    
    # Collect all atom indices that belong to the steroid nucleus.
    steroid_atoms = set()
    for idx in steroid_component:
        steroid_atoms.update(filtered_rings[idx])
    
    # As an extra precaution, check that the nucleus makes up at least 50% of the heavy atoms.
    if len(steroid_atoms) / mol.GetNumAtoms() < 0.5:
        return False, "Steroid nucleus constitutes less than 50% of the molecule's atoms"
    
    # Define SMARTS for a beta-oriented hydroxyl group.
    # In SMILES, a β OH group is typically expressed as [C@@H](O).
    beta_hydroxy_query = Chem.MolFromSmarts("[C@@H](O)")
    beta_matches = mol.GetSubstructMatches(beta_hydroxy_query)
    if not beta_matches:
        return False, "Beta-oriented hydroxyl group not found"
    
    # Check whether at least one beta–hydroxy carbon is part of the steroid nucleus.
    found_beta_in_nucleus = any(match[0] in steroid_atoms for match in beta_matches)
    if not found_beta_in_nucleus:
        return False, "Beta-oriented hydroxyl group not attached to the steroid nucleus"
    
    return True, "Molecule contains a steroid nucleus with a beta-oriented (3β) hydroxyl group"

# Example usage (uncomment to test):
# test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@@H](O)C(C)C"  # (24R)-24-hydroxycholesterol
# result, reason = is_3beta_hydroxy_steroid(test_smiles)
# print(result, reason)