"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5α-steroid
Definition: A 3-oxo steroid that has alpha configuration at position 5.
Heuristic improvements:
  • Filter rings to retain only 5- and 6-membered rings that are mostly carbon.
  • Among fused ring systems, search for exactly 4 rings with three 6-membered rings and one 5-membered ring.
  • Check that the union of the candidate nucleus atoms is in the typical range (15–20 atoms) and mostly carbon.
  • Confirm a ketone (C=O) is present on one of these nucleus atoms.
  • Require at least one chiral center in the nucleus with the '@@' (here, CHI_TETRAHEDRAL_CCW) configuration.
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
        bool: True if the molecule is classified as a 3-oxo-5α-steroid.
        str: Reason explaining the decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []
    for ring in ring_info:
        # Consider only rings with 5 or 6 atoms.
        if len(ring) in (5, 6):
            # Check that the ring is mostly carbon.
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            nC = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            if nC >= 0.8 * len(ring):
                candidate_rings.append((set(ring), len(ring)))
    if not candidate_rings:
        return False, "No suitable rings found; not steroid-like"

    # Build connectivity graph among these candidate rings.
    # Two rings are fused if they share at least 2 atoms.
    n = len(candidate_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(candidate_rings[i][0] & candidate_rings[j][0]) >= 2:
                graph[i].add(j)
                graph[j].add(i)

    # Find all connected components (sets of fused rings).
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

    # Try to find a valid 4-ring subset (steroid nucleus) in any component.
    steroid_component = None
    for comp in components:
        if len(comp) < 4:
            continue
        comp_list = list(comp)
        # If the entire connected component has exactly 4 rings, check it directly.
        if len(comp_list) == 4:
            sizes = [candidate_rings[i][1] for i in comp_list]
            if sizes.count(6) == 3 and sizes.count(5) == 1:
                # Obtain the union of atoms.
                nucleus_atoms = set()
                for i in comp_list:
                    nucleus_atoms |= candidate_rings[i][0]
                # Verify the nucleus has a realistic atom count and is mostly carbon.
                if 15 <= len(nucleus_atoms) <= 20:
                    nC_nucleus = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                    if nC_nucleus >= 0.7 * len(nucleus_atoms):
                        steroid_component = comp_list
                        break
        else:
            # More than 4 rings: search for a connected 4-ring combination.
            for subset in combinations(comp_list, 4):
                sizes = [candidate_rings[i][1] for i in subset]
                if not (sizes.count(6) == 3 and sizes.count(5) == 1):
                    continue
                # Check connectivity within the subset.
                subgraph = {i: graph[i] & set(subset) for i in subset}
                sub_visited = set()
                stack = [subset[0]]
                while stack:
                    cur = stack.pop()
                    if cur in sub_visited:
                        continue
                    sub_visited.add(cur)
                    stack.extend(subgraph[cur] - sub_visited)
                if len(sub_visited) == 4:
                    # Collect the union of atoms.
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
        return False, "Steroid nucleus not detected; fused set of 3 six-membered and 1 five-membered rings not found or not steroid-like"

    # Aggregate the atom indices in the detected steroid nucleus.
    nucleus_atoms = set()
    for idx in steroid_component:
        nucleus_atoms |= candidate_rings[idx][0]

    # Check for a ketone group (C=O) on an atom in the nucleus.
    ketone_pattern = Chem.MolFromSmarts("[#6;R]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_in_nucleus = False
    for match in ketone_matches:
        # match[0] is the carbonyl carbon.
        if match[0] in nucleus_atoms:
            atom = mol.GetAtomWithIdx(match[0])
            ring_neighbors = sum(1 for nb in atom.GetNeighbors() if nb.GetIdx() in nucleus_atoms)
            if ring_neighbors >= 2:
                ketone_in_nucleus = True
                break
    if not ketone_in_nucleus:
        return False, "Ketone group (C=O) not found in the steroid nucleus"

    # Check for an alpha configuration.
    # We expect at least one tetrahedral chiral center in the nucleus with the '@@' descriptor.
    alpha_found = False
    for idx in nucleus_atoms:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider tetrahedral centers.
        # In RDKit, chiral tags can be CHI_TETRAHEDRAL_CCW (often representing '@@') or CHI_TETRAHEDRAL_CW.
        if atom.GetChiralTag() in (Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, 
                                     Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW):
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                alpha_found = True
                break
    if not alpha_found:
        return False, "No chiral center with alpha (@@) configuration detected in the steroid nucleus"

    return True, "Molecule classified as a 3-oxo-5α-steroid (steroid nucleus, ketone, and alpha configuration detected)"