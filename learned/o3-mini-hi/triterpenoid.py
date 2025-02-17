"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid (any terpenoid derived from a triterpene)
Heuristic: Find the largest fused ring cluster and then “expand” its atom set by
including neighboring carbon atoms that are likely part of the core if bonded to
at least two atoms in the ring set. Then require that the (expanded) core contains
roughly 27–35 carbons, that there are at least 4 fused rings, and that the overall
molecular weight is consistent with a triterpenoid.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string using a heuristic.
    
    Steps:
      1. Parse the SMILES and extract all rings.
      2. Build a graph where each ring is a node (edge if rings share atoms) and find
         the largest connected (fused) ring cluster.
      3. The initial core is the union of atoms in the rings. Then, we expand the core by
         including any carbon atom outside the core if it touches at least two atoms already
         in the core. (This helps capture bridging carbons and parts of a rearranged skeleton.)
      4. Calculate the number of rings in the cluster and the number of carbon atoms in the
         expanded core.
      5. Apply the following criteria:
            - At least 4 fused rings.
            - The expanded core should have between 27 and 35 carbons.
            - The overall molecular weight should be above 400 Da.
            
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a triterpenoid, False otherwise.
        str: A message explaining the outcome.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check overall molecular weight: triterpenoids are generally > 400 Da.
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for triterpenoid"

    # Get all rings in the molecule as tuples of atom indices
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings detected in the structure"

    n_rings = len(all_rings)
    # Build connectivity among rings: rings are considered adjacent if they share at least one atom.
    adjacency = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if set(all_rings[i]).intersection(all_rings[j]):
                adjacency[i].add(j)
                adjacency[j].add(i)

    # Find connected components (fused clusters) of rings using depth-first search.
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

    # Choose the clustered set with the largest union of ring atoms.
    best_union = set()
    best_component = None
    for comp in components:
        union_atoms = set()
        for r in comp:
            union_atoms.update(all_rings[r])
        if len(union_atoms) > len(best_union):
            best_union = union_atoms
            best_component = comp

    if best_component is None:
        return False, "No fused ring cluster detected."

    # Count the number of fused rings in this cluster.
    n_fused_rings = len(best_component)

    # START: Expand the core.
    # We start with the union of all atoms in the fused ring cluster.
    core = set(best_union)
    added = True
    while added:
        added = False
        # Check all atoms not yet in the core.
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx in core:
                continue
            # Only consider carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Count how many neighbors are already in the core.
            nbrs = atom.GetNeighbors()
            count = sum(1 for nbr in nbrs if nbr.GetIdx() in core)
            # If at least two neighbors belong to the core, add this atom.
            if count >= 2:
                core.add(idx)
                added = True
    # END: Core expansion

    # Count carbons in the (expanded) core.
    n_carbons_in_core = sum(1 for idx in core if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)

    # Use a relaxed threshold to allow for some rearrangement.
    if n_fused_rings < 4:
        return False, f"Fused ring cluster has only {n_fused_rings} rings (expected at least 4 for a triterpenoid core)"
    if not (27 <= n_carbons_in_core <= 35):
        return False, f"Expanded fused ring core has {n_carbons_in_core} carbon atoms (expected between 27 and 35 for a triterpenoid)"
    
    return True, f"Fused ring core with {n_fused_rings} rings and {n_carbons_in_core} carbons is consistent with a triterpenoid"

# For testing purposes (this section may be removed or commented out in production)
if __name__ == "__main__":
    # Example: Deacetylnomilin (a triterpenoid) SMILES string
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)