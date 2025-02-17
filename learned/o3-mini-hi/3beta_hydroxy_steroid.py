"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 3β-hydroxy steroid
Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β–position.

This implementation parses the SMILES and then:
  1. Obtains rings (only 5- and 6-membered ones are candidates).
  2. Builds a graph of rings fused by at least 2 common atoms.
  3. Searches for a fused component (the steroid nucleus) that ideally:
       - Contains exactly 4 rings (the typical tetracyclic steroid core) where one is 5-membered and three are 6-membered;
         or if not available, then contains >=4 rings and the union of ring atoms is between 15 and 25 with >=70% carbons.
  4. Finally, it checks for at least one beta–oriented hydroxyl group ([C@@H](O))
     that is attached on a nucleus atom.
     
If all these conditions are met, the molecule is classified as a 3β–hydroxy steroid.
"""

from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3β–hydroxy steroid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Filter for rings of size 5 or 6.
    candidate_rings = []
    for ring in ring_info:
        if len(ring) in (5, 6):
            candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No 5- or 6-membered rings found"
    
    # Build a graph of rings fused by at least 2 atoms.
    ring_graph = {i:set() for i in range(len(candidate_rings))}
    for i in range(len(candidate_rings)):
        ring_i = set(candidate_rings[i])
        for j in range(i+1, len(candidate_rings)):
            ring_j = set(candidate_rings[j])
            if len(ring_i.intersection(ring_j)) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Discover connected components in the fuse graph.
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
            for neigh in ring_graph[node]:
                if neigh not in comp:
                    stack.append(neigh)
        seen.update(comp)
        fused_components.append(comp)
    
    steroid_nucleus = None
    nucleus_reason = ""
    # First try to find a fused component with exactly 4 rings and with 1 five-membered and 3 six-membered rings.
    for comp in fused_components:
        if len(comp) == 4:
            n5 = sum(1 for idx in comp if len(candidate_rings[idx]) == 5)
            n6 = sum(1 for idx in comp if len(candidate_rings[idx]) == 6)
            if n5 >= 1 and n6 >= 3:
                # Gather all atoms participating in these rings.
                nucleus_atoms = set()
                for idx in comp:
                    nucleus_atoms.update(candidate_rings[idx])
                steroid_nucleus = nucleus_atoms
                nucleus_reason = "Found tetracyclic (4-ring) nucleus with 1 five-membered and 3 six-membered rings"
                break

    # If no ideal tetracyclic nucleus was found, try candidates with >=4 rings and a nucleus size in [15,25],
    # and at least 70% carbons.
    if steroid_nucleus is None:
        for comp in fused_components:
            if len(comp) < 4:
                continue
            nucleus_atoms = set()
            for idx in comp:
                nucleus_atoms.update(candidate_rings[idx])
            if not (15 <= len(nucleus_atoms) <= 25):
                continue
            n_carbons = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if (n_carbons / len(nucleus_atoms)) < 0.70:
                continue
            # Additionally, require at least one 5-membered and one 6-membered ring overall.
            n5 = sum(1 for idx in comp if len(candidate_rings[idx]) == 5)
            n6 = sum(1 for idx in comp if len(candidate_rings[idx]) == 6)
            if n5 < 1 or n6 < 1:
                continue
            steroid_nucleus = nucleus_atoms
            nucleus_reason = f"Found fused nucleus (from {len(comp)} rings, union size {len(nucleus_atoms)}, carbon fraction {n_carbons/len(nucleus_atoms):.2f})"
            break

    if steroid_nucleus is None:
        return False, "Steroid nucleus not found (no appropriate fused tetracyclic cluster detected)"
    
    # (Optional) Check that the nucleus is a reasonable fraction of the whole molecule.
    if len(steroid_nucleus)/mol.GetNumAtoms() < 0.20:
        return False, "Nucleus is too small relative to the molecule"
    
    # Look for beta–oriented hydroxyl group.
    # SMARTS [C@@H](O) looks for a chiral carbon with an –OH.
    beta_oh_query = Chem.MolFromSmarts("[C@@H](O)")
    beta_matches = mol.GetSubstructMatches(beta_oh_query)
    if not beta_matches:
        return False, "No beta–oriented hydroxyl group ([C@@H](O)) found"
    
    # Require that at least one beta-OH is attached to an atom in the nucleus.
    for match in beta_matches:
        # match[0] gives the index of the chiral carbon (with –OH attached).
        if match[0] in steroid_nucleus:
            return True, "Molecule contains a steroid nucleus (" + nucleus_reason + ") with a beta-oriented (3β) hydroxyl group"
    
    return False, "No beta–oriented hydroxyl group found on the steroid nucleus"

# Example usage:
# test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@@H](O)C(C)C"
# result, reason = is_3beta_hydroxy_steroid(test_smiles)
# print(result, reason)