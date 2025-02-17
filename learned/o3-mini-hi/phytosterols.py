"""
Classifies: CHEBI:26125 phytosterols
"""
#!/usr/bin/env python
"""
Classifies: Phytosterols – Sterols similar to cholesterol that occur in plants.
Phytosterols contain a characteristic steroid nucleus (a fused tetracyclic system with three 6-membered rings and one 5-membered ring),
often a hydroxyl group, and generally have molecular weights in the 300-600 Da range.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES string.
      2. Extracts all ring information.
      3. Builds a connectivity graph of rings (rings that are fused share at least 2 atoms).
      4. Searches for any connected component of rings that has at least 4 fused rings,
         with at least one 5-membered ring and at least three 6-membered rings.
         (This is a heuristic for the steroid nucleus.)
      5. Checks for at least one hydroxyl group ([OX2H]) in the molecule.
      6. Verifies the molecular weight is within the 300-600 Da range.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phytosterol, False otherwise.
        str: Explanation (reason) for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of tuples containing atom indices for each ring
    if not rings:
        return False, "No rings found; not a steroid structure"

    # Build a graph of rings: nodes are ring indices; connect two rings if they share at least 2 atoms (fused)
    n_rings = len(rings)
    ring_graph = {i: [] for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(rings[i])
        for j in range(i+1, n_rings):
            set_j = set(rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_graph[i].append(j)
                ring_graph[j].append(i)

    # Helper function: get connected component starting at a given node (ring index)
    def dfs(start, visited):
        stack = [start]
        comp = set()
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                comp.add(node)
                for neigh in ring_graph[node]:
                    if neigh not in visited:
                        stack.append(neigh)
        return comp

    visited = set()
    steroid_found = False
    reason_detail = ""
    # Look for any connected component of rings that may represent a steroid nucleus
    for i in range(n_rings):
        if i in visited:
            continue
        comp = dfs(i, visited)
        if len(comp) < 4:
            # not enough fused rings
            continue
        # For each ring in the component, count how many are 5-membered and how many are 6-membered
        count5 = 0
        count6 = 0
        # Also count carbon atoms in these rings (avoiding duplicates) as a heuristic check
        comp_atom_indices = set()
        for idx in comp:
            comp_atom_indices.update(rings[idx])
            ring_size = len(rings[idx])
            if ring_size == 5:
                count5 += 1
            elif ring_size == 6:
                count6 += 1
        # In a classic steroid nucleus we expect at least 1 five-membered and 3 six-membered rings.
        if count5 >= 1 and count6 >= 3:
            # Additionally, the steroid nucleus should contain a significant number of carbons.
            carbon_count = sum(1 for a in mol.GetAtoms() if a.GetIdx() in comp_atom_indices and a.GetAtomicNum() == 6)
            if carbon_count >= 15:  # expecting around 17 carbons in the nucleus; allow a little flexibility
                steroid_found = True
                break

    if not steroid_found:
        return False, "Molecule does not contain a characteristic steroid nucleus (fused 3 six-membered and 1 five-membered rings)"

    # Check for the presence of at least one hydroxyl group (–OH), common in phytosterols.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxyl_pattern is None:
        return False, "Error creating hydroxyl SMARTS pattern"
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Molecule lacks a hydroxyl group typical of phytosterols"

    # Check for molecular weight (typical sterol range ca. 300-600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the typical sterol range for phytosterols"

    return True, "Molecule exhibits a fused steroid nucleus, a hydroxyl group, and appropriate molecular weight, consistent with a phytosterol"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Example: ergosta-5,7-dien-3beta-ol
    test_smiles = "[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@H](C)C(C)C"
    result, reason = is_phytosterols(test_smiles)
    print(result, reason)