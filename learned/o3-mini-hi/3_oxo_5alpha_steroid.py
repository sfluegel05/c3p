"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5α-steroid
Definition: A 3-oxo steroid that has alpha configuration at position 5.
Heuristic: 
  • The molecule must have a fused ring system containing exactly 4 rings that are connected via sharing at least 2 atoms.
    In a typical steroid nucleus there are three six-membered rings and one five-membered ring.
  • At least one ketone (C=O) group must be present and its carbon must lie in one of these rings.
  • At least one chiral center in the fused set is marked with the '@@' stereochemistry (heuristically taken as alpha).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5α-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-oxo-5α-steroid.
        str: Reason explaining the result.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information for all rings.
    ring_info = mol.GetRingInfo().AtomRings()
    # We are interested only in rings of size 5 or 6.
    candidate_rings = []
    for ring in ring_info:
        if len(ring) in (5, 6):
            candidate_rings.append( (set(ring), len(ring)) )
    
    if not candidate_rings:
        return False, "No rings of size 5 or 6 found; cannot be steroid-like"
    
    # Build a connectivity graph amongst these rings.
    # Two rings are considered fused if they share at least 2 atoms.
    n = len(candidate_rings)
    # Graph: index -> set of connected ring indices.
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(candidate_rings[i][0] & candidate_rings[j][0]) >= 2:
                graph[i].add(j)
                graph[j].add(i)
    
    # Find connected components (groups of fused rings).
    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                cur = stack.pop()
                if cur in comp:
                    continue
                comp.add(cur)
                for neigh in graph[cur]:
                    if neigh not in comp:
                        stack.append(neigh)
            visited |= comp
            components.append(comp)
    
    # Among all connected components, check if any has exactly 4 rings with 3 six-membered and 1 five-membered rings.
    steroid_component = None
    for comp in components:
        if len(comp) < 4:
            continue  # too few rings in this fused set
        # In some molecules extra rings may be present; so we search for any subset of 4 rings that are mutually connected.
        # For simplicity, if the connected component has exactly 4 rings, use it.
        if len(comp) == 4:
            sizes = [candidate_rings[i][1] for i in comp]
            if sizes.count(6) == 3 and sizes.count(5) == 1:
                steroid_component = comp
                break
        # If more than 4 rings are present, try to find a valid 4-ring combination.
        else:
            comp_list = list(comp)
            from itertools import combinations
            for subset in combinations(comp_list, 4):
                sizes = [candidate_rings[i][1] for i in subset]
                # To ensure they are mutually fused, check that for every pair in the subset, the rings are fused (share >=2 atoms) or are connected through others.
                # Here we do a simple check on the induced subgraph connectivity.
                subgraph = {i: graph[i] & set(subset) for i in subset}
                # Do a DFS in the subgraph
                sub_visited = set()
                stack = [subset[0]]
                while stack:
                    cur = stack.pop()
                    if cur in sub_visited:
                        continue
                    sub_visited.add(cur)
                    stack.extend(subgraph[cur] - sub_visited)
                if len(sub_visited) == 4 and sizes.count(6) == 3 and sizes.count(5) == 1:
                    steroid_component = set(subset)
                    break
            if steroid_component is not None:
                break

    if steroid_component is None:
        return False, "Steroid nucleus not detected; fused set of 3 six-membered and 1 five-membered rings not found"
    
    # Get the union of atom indices that belong to the steroid nucleus.
    steroid_atoms = set()
    for i in steroid_component:
        steroid_atoms |= candidate_rings[i][0]
    
    # Check for the 3-oxo group, a ketone (C=O) that is part of a ring in the steroid nucleus.
    ketone_pattern = Chem.MolFromSmarts("[#6;R]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_in_nucleus = False
    for match in ketone_matches:
        # match[0] is the carbonyl carbon.
        if match[0] in steroid_atoms:
            # Further requirement: the carbonyl should be attached to at least 2 ring atoms from the nucleus.
            atom = mol.GetAtomWithIdx(match[0])
            ring_neighbors = sum(1 for nb in atom.GetNeighbors() if nb.GetIdx() in steroid_atoms)
            if ring_neighbors >= 2:
                ketone_in_nucleus = True
                break
    if not ketone_in_nucleus:
        return False, "Ketone group (C=O) not found in the steroid nucleus"

    # Check for the alpha configuration.
    # We expect at least one chiral center in the nucleus marked with the '@@' descriptor.
    alpha_found = False
    for idx in steroid_atoms:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider tetrahedral chiral centers.
        if atom.GetChiralTag() in (Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW):
            # Heuristically choose CCW as the alpha configuration.
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                alpha_found = True
                break
    if not alpha_found:
        return False, "No chiral center with alpha (@@) configuration detected in the steroid nucleus"

    return True, "Molecule classified as a 3-oxo-5α-steroid (fused steroid nucleus with ketone and alpha configuration detected)"
    
# The function can be tested by calling is_3_oxo_5alpha_steroid(smiles) with example SMILES.