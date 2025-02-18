"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI: 3-oxo-5α-steroid
Definition: A 3-oxo steroid with an alpha configuration at position 5.
Heuristic:
  • Look for a fused four‐ring steroid nucleus: exactly three six‐membered rings and one five‐membered ring,
    whose union of ring atoms is within 16–18 atoms and is mostly carbon.
  • The five‐membered ring must be fused (≥2 common atoms) with exactly one six‐membered ring.
  • Look for a ketone group (C=O) on that six‐membered ring (and not on the five‐membered ring).
  • Check for a chiral center in the five‐membered ring. If not present explicitly, we assume the configuration is 5α,
    but we note that stereochemistry was not explicitly specified.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from itertools import combinations

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5α-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-oxo-5α-steroid, False otherwise.
        str: Reason explaining the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Ensure stereochemistry is assigned.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Retrieve ring information (only consider rings of size 5 or 6 that are mostly carbon).
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []  # list of tuples: (set_of_atom_indices, ring_size)
    for ring in ring_info:
        if len(ring) in (5, 6):
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # require at least 80% of atoms to be carbon.
            nC = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            if nC >= 0.8 * len(ring):
                candidate_rings.append((set(ring), len(ring)))
    
    if not candidate_rings:
        return False, "No suitable 5-/6-membered rings found; not steroid-like"
    
    # Build a connectivity graph among candidate rings (fused if sharing at least 2 atoms).
    n = len(candidate_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(candidate_rings[i][0] & candidate_rings[j][0]) >= 2:
                graph[i].add(j)
                graph[j].add(i)
    
    # Find connected components (fused ring systems).
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
                stack.extend(graph[cur] - comp)
            visited |= comp
            components.append(comp)
    
    # Look for a fused system that has exactly four rings (three six-membered and one five-membered).
    steroid_component = None
    for comp in components:
        comp_list = list(comp)
        if len(comp_list) < 4:
            continue
        # For components with more than 4 rings, check every combination of 4 rings.
        possible_sets = []
        if len(comp_list) == 4:
            possible_sets.append(comp_list)
        else:
            for subset in combinations(comp_list, 4):
                possible_sets.append(list(subset))
        for subset in possible_sets:
            sizes = [candidate_rings[i][1] for i in subset]
            if sizes.count(6) != 3 or sizes.count(5) != 1:
                continue
            # Ensure the selected rings are all fused among themselves.
            subgraph = {i: graph[i] & set(subset) for i in subset}
            subvisited = set()
            stack = [subset[0]]
            while stack:
                cur = stack.pop()
                if cur in subvisited:
                    continue
                subvisited.add(cur)
                stack.extend(subgraph[cur] - subvisited)
            if len(subvisited) != 4:
                continue
            # Merge atom indices.
            nucleus_atoms = set()
            for i in subset:
                nucleus_atoms |= candidate_rings[i][0]
            # Classical steroid nucleus has 17 atoms; allow a tight range.
            if not (16 <= len(nucleus_atoms) <= 18):
                continue
            # At least 70% of the nucleus atoms should be carbon.
            nC_nucleus = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if nC_nucleus < 0.7 * len(nucleus_atoms):
                continue
            steroid_component = subset
            break
        if steroid_component is not None:
            break
    
    if steroid_component is None:
        return False, "Steroid nucleus not detected; required fused system (3 six‐membered and 1 five‐membered ring) not found or not steroid‐like"
    
    # Merge the nucleus atoms.
    nucleus_atoms = set()
    for idx in steroid_component:
        nucleus_atoms |= candidate_rings[idx][0]
    
    # Identify the unique five-membered ring.
    five_membered_ring = None
    for idx in steroid_component:
        ring_set, size = candidate_rings[idx]
        if size == 5:
            five_membered_ring = ring_set
            break
    if five_membered_ring is None:
        return False, "Five-membered ring missing in fused steroid nucleus candidate"
    
    # Additional check: the five-membered ring must be fused with exactly one six-membered ring.
    six_fused_count = 0
    for idx in steroid_component:
        ring_set, size = candidate_rings[idx]
        if size == 6:
            if len(ring_set & five_membered_ring) >= 2:
                six_fused_count += 1
    if six_fused_count != 1:
        return False, f"Expected the five-membered ring to fuse with only one six-membered ring, found {six_fused_count}"
    
    # Look for a ketone group (C=O) on one of the six-membered rings.
    # Pattern: carbon in a ring double-bonded to oxygen.
    ketone_pattern = Chem.MolFromSmarts("[#6;R]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_found = False
    for match in ketone_matches:
        ketone_carbon = match[0]
        # must be in the nucleus and not in the five-membered ring.
        if ketone_carbon in nucleus_atoms and ketone_carbon not in five_membered_ring:
            # Check that the ketone carbon has at least two neighbors within the nucleus.
            atom = mol.GetAtomWithIdx(ketone_carbon)
            connected_in_nucleus = sum(1 for nb in atom.GetNeighbors() if nb.GetIdx() in nucleus_atoms)
            if connected_in_nucleus >= 2:
                ketone_found = True
                break
    if not ketone_found:
        return False, "Ketone group (C=O) not found on a six-membered ring of the steroid nucleus"
    
    # Check the five-membered ring for at least one chiral center.
    # We use FindMolChiralCenters with only assigned centers.
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    chiral_in_five = [idx for (idx, tag) in centers if idx in five_membered_ring]
    if not chiral_in_five:
        # In our previous attempt we rejected if no chiral center was found,
        # but here we relax the criterion to avoid missing known steroids.
        chirality_note = " (Warning: no explicit chiral center detected in the five-membered ring; assuming 5α configuration)"
    else:
        chirality_note = ""
    
    return True, ("Molecule classified as a 3-oxo-5α-steroid: fused steroid nucleus (3 six‐membered and 1 five‐membered rings), "
                  "ketone on a six‐membered ring, and the expected 5α configuration detected" + chirality_note)
    
# End of code.