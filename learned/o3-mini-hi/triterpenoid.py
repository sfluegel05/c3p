"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid (any terpenoid derived from a triterpene)
The heuristic requires a fused polycyclic core (typically 4–6 rings) whose union contains about 30 carbons (we allow 27–33).
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    
    The method uses a heuristic:
      1. The molecule is parsed and its ring system is extracted.
      2. The rings are grouped into fused clusters (rings sharing atoms).
      3. The largest fused ring cluster is assumed to be the triterpenoid core.
      4. If the union of atoms in this cluster contains about 27–33 carbons
         and the cluster comprises at least 4 rings, we classify the molecule as a triterpenoid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a triterpenoid, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings in the molecule as tuples of atom indices
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings detected in the structure"

    # Group rings into clusters that are fused (i.e. share at least one atom)
    # We treat each ring as a node in a graph and add an edge if two rings share at least one atom.
    n_rings = len(all_rings)
    # Build an adjacency list for ring indices
    adjacency = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if set(all_rings[i]).intersection(all_rings[j]):
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Now, find connected components in the ring graph using DFS.
    visited = set()
    components = []
    for i in range(n_rings):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                comp.add(node)
                stack.extend(adjacency[node] - visited)
            components.append(comp)
    
    # For each component, take the union of atom indices and count carbons.
    def count_carbon_atoms(atom_indices):
        count = 0
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                count += 1
        return count

    # Choose the component that has the largest union of ring atoms.
    best_union = set()
    best_component = None
    for comp in components:
        union_atoms = set()
        for r in comp:
            union_atoms.update(all_rings[r])
        if len(union_atoms) > len(best_union):
            best_union = union_atoms
            best_component = comp

    # Count the number of rings in this fused system and carbons in the union.
    n_fused_rings = len(best_component) if best_component is not None else 0
    n_carbons_in_core = count_carbon_atoms(best_union)

    # Heuristic: a typical triterpenoid core is derived from a C30 framework,
    # sometimes featuring a slightly reduced amount (27–33 carbons) and 4 or more fused rings.
    if n_fused_rings < 4:
        return False, f"Fused ring cluster has only {n_fused_rings} rings (expected at least 4 for a triterpenoid core)"
    if not (27 <= n_carbons_in_core <= 33):
        return False, f"Fused ring core has {n_carbons_in_core} carbon atoms (expected between 27 and 33 for a triterpenoid)"
    
    return True, f"Fused ring core with {n_fused_rings} rings and {n_carbons_in_core} carbons is consistent with a triterpenoid"

# For testing purposes (feel free to remove or comment out in production)
if __name__ == "__main__":
    # Example: Deacetylnomilin (a triterpenoid) SMILES string
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)