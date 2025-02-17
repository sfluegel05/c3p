"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5α-steroid
Definition: A 3-oxo steroid that has alpha configuration at position 5.
Heuristic improvements in this version:
  • Only consider fused sets of exactly 4 rings (3 six-membered and 1 five-membered)
  • Check that the union of rings (the steroid nucleus) contains 15–20 atoms and is mostly carbon.
  • Confirm that a ketone (C=O) is present on a nucleus atom.
  • Require that the unique 5-membered ring has at least one chiral center with the CHI_TETRAHEDRAL_CCW ("@@") configuration.
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

    # Get ring information and filter for rings of size 5 or 6 that are mostly carbon.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []  # list of tuples: (set_of_atom_indices, ring_size)
    for ring in ring_info:
        if len(ring) in (5, 6):
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            nC = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            if nC >= 0.8 * len(ring):
                candidate_rings.append((set(ring), len(ring)))
    if not candidate_rings:
        return False, "No suitable rings found; not steroid-like"

    # Build connectivity graph among candidate rings. Two rings are fused if they share at least 2 atoms.
    n = len(candidate_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(candidate_rings[i][0] & candidate_rings[j][0]) >= 2:
                graph[i].add(j)
                graph[j].add(i)

    # Find connected components (i.e. fused ring systems)
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

    # Look for a valid four-ring fused system that meets steroid nucleus criteria.
    steroid_component = None
    for comp in components:
        comp_list = list(comp)
        if len(comp_list) < 4:
            continue
        # If the entire component has exactly four rings, check it directly.
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
            # If more than four rings are fused, try every combination of 4.
            for subset in combinations(comp_list, 4):
                sizes = [candidate_rings[i][1] for i in subset]
                if sizes.count(6) != 3 or sizes.count(5) != 1:
                    continue
                # Check that the selected rings are themselves fused.
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
        return False, "Steroid nucleus not detected; required fused system (3 six‐membered and 1 five‐membered rings) not found or not steroid-like"

    # Combine the atoms from the selected rings (the steroid nucleus).
    nucleus_atoms = set()
    for idx in steroid_component:
        nucleus_atoms |= candidate_rings[idx][0]

    # Check for a ketone group (C=O) within the nucleus.
    ketone_pattern = Chem.MolFromSmarts("[#6;R]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_in_nucleus = False
    for match in ketone_matches:
        # Using the carbonyl carbon index.
        if match[0] in nucleus_atoms:
            # Ensure that this carbon is bonded to at least two atoms in the nucleus.
            atom = mol.GetAtomWithIdx(match[0])
            ring_nb = sum(1 for nb in atom.GetNeighbors() if nb.GetIdx() in nucleus_atoms)
            if ring_nb >= 2:
                ketone_in_nucleus = True
                break
    if not ketone_in_nucleus:
        return False, "Ketone group (C=O) not found in the steroid nucleus"

    # Identify the unique five-membered ring in the nucleus candidate.
    five_membered_ring = None
    for idx in steroid_component:
        ring_set, size = candidate_rings[idx]
        if size == 5:
            five_membered_ring = ring_set
            break
    if five_membered_ring is None:
        return False, "Five-membered ring missing in fused steroid nucleus candidate"

    # Check that at least one atom in the 5-membered ring bears the CHI_TETRAHEDRAL_CCW chiral tag ("@@").
    alpha_found = False
    for idx in five_membered_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            alpha_found = True
            break
    if not alpha_found:
        return False, "No chiral center with alpha (@@) configuration detected in the five-membered ring"

    return True, "Molecule classified as a 3-oxo-5α-steroid (steroid nucleus, ketone, and alpha configuration detected)"