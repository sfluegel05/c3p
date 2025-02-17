"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3β-hydroxy steroid
Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β‐position.
Our improved classifier proceeds as follows:
  1. We parse the SMILES and extract rings (only rings of size 5 or 6 that are mostly carbon).
  2. We build a “fused ring” graph (rings sharing at least 2 atoms) and look for a fused cluster 
     containing at least 4 rings – ideally one 5-membered ring and three 6-membered rings.
  3. We collect the atoms making up that candidate steroid nucleus.
  4. We require that the nucleus makes up at least a modest fraction (~20%) of the heavy atoms.
  5. We then search for beta–oriented hydroxyl groups (SMARTS “[C@@H](O)”) and require that at least one 
     such beta hydroxyl is attached to an atom that is part of the nucleus.
If all these conditions are met, we classify the molecule as a 3β–hydroxy steroid.
"""

from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.
    
    The classifier checks:
      1. Whether a steroid nucleus is present – we look for a fused tetracyclic system
         composed of rings of size 5 or 6. We filter out rings which are not predominantly 
         carbon (>=67% C atoms).
      2. That we can find a connected fused ring cluster containing at least four rings,
         with at least one 5-membered ring and three 6-membered rings.
      3. That the nucleus spans at least 20% of all heavy atoms (to avoid cases where the steroid 
         is only a tiny fragment of a hugely decorated molecule).
      4. That there is at least one beta–oriented hydroxyl group ([C@@H](O)) attached to an atom 
         in the nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a 3β-hydroxy steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Filter rings: we allow only rings of size 5 or 6 having at least 67% carbons.
    filtered_rings = []
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        n_carbon = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbon / len(ring) >= 0.67:
            filtered_rings.append(ring)
    
    if not filtered_rings:
        return False, "No rings of size 5 or 6 with sufficient carbon content found"
    
    # Build a graph of rings:
    # Two rings are considered fused if they share at least 2 atoms.
    ring_graph = {i: set() for i in range(len(filtered_rings))}
    for i in range(len(filtered_rings)):
        set_i = set(filtered_rings[i])
        for j in range(i+1, len(filtered_rings)):
            set_j = set(filtered_rings[j])
            if len(set_i & set_j) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components (clusters of fused rings).
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
    
    # Look for a fused component with at least 4 rings,
    # and ideally having at least one 5-membered ring and three 6-membered rings.
    steroid_comp = None
    for comp in fused_components:
        if len(comp) < 4:
            continue
        n5 = sum(1 for idx in comp if len(filtered_rings[idx]) == 5)
        n6 = sum(1 for idx in comp if len(filtered_rings[idx]) == 6)
        if n5 >= 1 and n6 >= 3:
            steroid_comp = comp
            break
    if steroid_comp is None:
        return False, "Steroid nucleus not found (no fused tetracyclic cluster with 1 five‐membered and 3 six‐membered rings)"
    
    # Get the set of all atoms that belong to the steroid nucleus.
    nucleus_atoms = set()
    for idx in steroid_comp:
        nucleus_atoms.update(filtered_rings[idx])
    
    # Optionally, require that the nucleus represents at least 20% of the heavy atoms in the molecule.
    total_atoms = mol.GetNumAtoms()
    if len(nucleus_atoms) / total_atoms < 0.20:
        return False, "Steroid nucleus constitutes less than 20% of molecule's atoms"
    
    # Search for beta-oriented hydroxyl groups.
    # The SMARTS for a beta-oriented hydroxyl group in our context is [C@@H](O)
    beta_oh_query = Chem.MolFromSmarts("[C@@H](O)")
    beta_matches = mol.GetSubstructMatches(beta_oh_query)
    if not beta_matches:
        return False, "No beta-oriented hydroxyl group ([C@@H](O)) found"
    
    # Check whether at least one beta –OH is attached to an atom that is part of the nucleus.
    for match in beta_matches:
        # match[0] is the carbon with defined stereochemistry bearing an -OH.
        if match[0] in nucleus_atoms:
            return True, "Molecule contains a steroid nucleus with a beta-oriented (3β) hydroxyl group"
    
    return False, "Beta-oriented hydroxyl group not attached to the steroid nucleus"


# Example usage (uncomment to test):
# test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@@H](O)C(C)C"  # (24R)-24-hydroxycholesterol
# result, reason = is_3beta_hydroxy_steroid(test_smiles)
# print(result, reason)