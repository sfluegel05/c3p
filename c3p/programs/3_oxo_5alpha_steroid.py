"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI: 3-oxo-5α-steroid
Definition: A 3-oxo steroid with an alpha configuration at position 5.
Heuristic:
  • Look for a fused four‐ring steroid nucleus: exactly three six‐membered rings and one five‐membered ring that are fused (sharing ≥2 atoms) and whose union has 15–20 atoms and is mostly carbon.
  • Require the presence of a ketone group (C=O) on a ring atom of the nucleus but not on the unique five‐membered ring.
  • Confirm that the unique five‐membered ring has at least one explicitly chiral center (regardless of CW or CCW) to indicate a defined configuration at position 5.
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
    
    # Retrieve ring information: consider rings of size 5 or 6 that are mostly carbon.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []  # list of tuples: (set_of_atom_indices, ring_size)
    for ring in ring_info:
        if len(ring) in (5, 6):
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Require at least 80% carbon atoms in the ring.
            nC = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            if nC >= 0.8 * len(ring):
                candidate_rings.append((set(ring), len(ring)))
    if not candidate_rings:
        return False, "No suitable rings found; not steroid-like"
    
    # Build a connectivity graph among candidate rings.
    # Two rings are fused if they share at least 2 atoms.
    n = len(candidate_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(candidate_rings[i][0] & candidate_rings[j][0]) >= 2:
                graph[i].add(j)
                graph[j].add(i)
    
    # Find connected components of fused rings.
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
    
    # Look for a fused system that has exactly four rings: three six-membered and one five-membered.
    steroid_component = None
    for comp in components:
        comp_list = list(comp)
        if len(comp_list) < 4:
            continue
        # Check either the entire component if it has exactly 4 rings...
        if len(comp_list) == 4:
            sizes = [candidate_rings[i][1] for i in comp_list]
            if sizes.count(6) == 3 and sizes.count(5) == 1:
                # Merge all atoms from these rings.
                nucleus_atoms = set()
                for i in comp_list:
                    nucleus_atoms |= candidate_rings[i][0]
                if 15 <= len(nucleus_atoms) <= 20:
                    nC_nucleus = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                    if nC_nucleus >= 0.7 * len(nucleus_atoms):
                        steroid_component = comp_list
                        break
        else:
            # If there are more fused rings then try every combination of 4.
            for subset in combinations(comp_list, 4):
                sizes = [candidate_rings[i][1] for i in subset]
                if sizes.count(6) != 3 or sizes.count(5) != 1:
                    continue
                # Ensure the selected four rings are all fused.
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
                # Check the union of atoms.
                nucleus_atoms = set()
                for i in subset:
                    nucleus_atoms |= candidate_rings[i][0]
                if not (15 <= len(nucleus_atoms) <= 20):
                    continue
                nC_nucleus = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                if nC_nucleus < 0.7 * len(nucleus_atoms):
                    continue
                steroid_component = list(subset)
                break
        if steroid_component is not None:
            break

    if steroid_component is None:
        return False, "Steroid nucleus not detected; required fused system (3 six‐membered and 1 five‐membered rings) not found or not steroid‐like"
    
    # Combine the nucleus atoms.
    nucleus_atoms = set()
    for idx in steroid_component:
        nucleus_atoms |= candidate_rings[idx][0]
    
    # Identify the unique five-membered ring in the nucleus.
    five_membered_ring = None
    for idx in steroid_component:
        ring_set, size = candidate_rings[idx]
        if size == 5:
            five_membered_ring = ring_set
            break
    if five_membered_ring is None:
        return False, "Five-membered ring missing in fused steroid nucleus candidate"
    
    # Look for a ketone group (C=O) on the nucleus. 
    # But require that the carbonyl carbon is NOT part of the five-membered ring, 
    # because the 3-oxo group should be on one of the six-membered rings (usually ring A).
    ketone_pattern = Chem.MolFromSmarts("[#6;R]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_in_six_ring = False
    for match in ketone_matches:
        ketone_idx = match[0]  # the carbon in C=O
        if ketone_idx in nucleus_atoms and ketone_idx not in five_membered_ring:
            # Check that the ketone carbon is well integrated: at least two bonds to the nucleus.
            atom = mol.GetAtomWithIdx(ketone_idx)
            connected_in_nucleus = sum(1 for nb in atom.GetNeighbors() if nb.GetIdx() in nucleus_atoms)
            if connected_in_nucleus >= 2:
                ketone_in_six_ring = True
                break
    if not ketone_in_six_ring:
        return False, "Ketone group (C=O) not found on a six-membered ring of the steroid nucleus"
    
    # Check the five-membered ring for at least one chiral center.
    # Instead of strictly demanding CHI_TETRAHEDRAL_CCW, we accept any explicit chiral tag.
    alpha_found = False
    for idx in five_membered_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            alpha_found = True
            break
    if not alpha_found:
        return False, "No chiral center detected in the five-membered ring (alpha configuration not established)"
    
    return True, "Molecule classified as a 3-oxo-5α-steroid (steroid nucleus with 3 six‐membered and 1 five‐membered rings, ketone on a six‐membered ring, and defined chiral center in the five‐membered ring detected)"