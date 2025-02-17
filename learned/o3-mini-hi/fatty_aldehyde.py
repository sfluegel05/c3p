"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehyde
Definition: An aldehyde formally arising from reduction of the carboxylic acid of its corresponding fatty acid,
with a long acyclic chain that is mainly carbon/hydrogen and a terminal aldehyde group.

Improvements over the previous version:
  • Only allows elements C, H, and O.
  • Requires an acyclic molecular structure.
  • Demands at least 6 carbon atoms.
  • Requires that the non‐hydrogen atoms are predominantly carbon (fraction ≥ 0.75).
  • Requires exactly one overall carbonyl group ([CX3]=[OX1]) in the molecule.
  • The aldehyde carbon in that carbonyl must be terminal (attached to exactly one carbon and not in a ring).
  • The longest carbon chain (by connectivity among carbons) must encompass at least 75% of all C atoms;
    furthermore, the terminal aldehyde must lie at one end of that chain.
  • The chain must not be excessively conjugated – measured here as having a low fraction of sp² carbons (excluding the carbonyl)
    in the chain (threshold set to 0.3).
"""

from rdkit import Chem
from collections import deque
from rdkit.Chem import rdchem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule qualifies as a fatty aldehyde.
    
    The method proceeds as follows:
      1. Parse the SMILES string.
      2. Confirm that only allowed elements (C, H, O) are present.
      3. Ensure that the molecule is acyclic.
      4. Verify that there are at least 6 carbon atoms.
      5. Confirm that heavy atoms (non-H) are predominantly carbons (fraction ≥ 0.75).
      6. Ensure that the molecule contains exactly one carbonyl group ([CX3]=[OX1]).
      7. Confirm that the carbonyl appears as a terminal aldehyde:
         - The aldehyde carbon must have one bonded carbon (and at least one hydrogen via [CX3H1](=O))
         - It must not be in a ring.
      8. Calculate the connected carbon subgraph and determine the longest chain.
         - The longest chain must account for at least 75% of carbon atoms.
         - The terminal aldehyde carbon must be at one end of that chain.
      9. Finally, assess the saturation of the fatty chain.
         - Excluding the aldehyde carbon, the fraction of sp²-hybridized carbons must be below 0.3.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple where the boolean indicates acceptance and the string explains the decision.
    """
    # 1. Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check that only allowed elements are present: C (6), H (1), and O (8)
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: {atom.GetSymbol()}"
    
    # 3. Ensure that the molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), which is not typical for a fatty aldehyde"
    
    # 4. Count carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbon_atoms)
    if n_carbons < 6:
        return False, f"Not enough carbon atoms ({n_carbons} found; need at least 6) for a fatty chain"
    
    # 5. Check heavy atom composition (non-H atoms).
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return False, "No heavy atoms found"
    carbon_fraction = n_carbons / len(heavy_atoms)
    if carbon_fraction < 0.75:
        return False, f"Molecule has a low carbon fraction ({carbon_fraction:.2f}); might have additional functionalities"
    
    # 6. Identify the carbonyl group: exactly one [CX3]=[OX1]
    carbonyl_smarts = "[CX3]=[OX1]"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) != 1:
        return False, f"Expected exactly one carbonyl group; found {len(carbonyl_matches)}"
    
    # 7. Confirm terminal aldehyde functionality with at least one hydrogen on the carbonyl C.
    aldehyde_smarts = "[CX3H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    terminal_aldehyde_idx = None
    terminal_aldehyde_count = 0
    for match in aldehyde_matches:
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        # Exclude if the carbonyl C is in a ring.
        if aldehyde_c.IsInRing():
            continue
        # Count neighboring carbon atoms (ignore oxygen neighbors).
        carbon_neighbors = [nbr for nbr in aldehyde_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:  # terminal condition
            terminal_aldehyde_count += 1
            terminal_aldehyde_idx = match[0]
    if terminal_aldehyde_count != 1:
        return False, f"Expected exactly one terminal aldehyde group; found {terminal_aldehyde_count}"
    
    # 8. Build the carbon connectivity graph.
    carbon_indices = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    carbon_graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        if idx1 in carbon_indices and idx2 in carbon_indices:
            carbon_graph[idx1].append(idx2)
            carbon_graph[idx2].append(idx1)
    
    # Confirm that all carbon atoms form a single connected component.
    visited = set()
    def bfs(start):
        q = deque([start])
        comp = set()
        while q:
            cur = q.popleft()
            if cur in comp:
                continue
            comp.add(cur)
            for nei in carbon_graph[cur]:
                if nei not in comp:
                    q.append(nei)
        return comp
    comp = bfs(next(iter(carbon_indices)))
    if len(comp) != len(carbon_indices):
        return False, "Carbon atoms are not all connected in one chain-like structure"
    
    # Calculate the tree diameter (longest chain) in the acyclic carbon graph.
    def bfs_farthest(start):
        seen = {start: 0}
        q = deque([start])
        parent = {start: None}
        farthest_node = start
        max_dist = 0
        while q:
            cur = q.popleft()
            for nei in carbon_graph[cur]:
                if nei not in seen:
                    seen[nei] = seen[cur] + 1
                    parent[nei] = cur
                    if seen[nei] > max_dist:
                        max_dist = seen[nei]
                        farthest_node = nei
                    q.append(nei)
        return farthest_node, max_dist, parent

    # First BFS: from an arbitrary carbon.
    arbitrary = next(iter(carbon_indices))
    far_node, _, _ = bfs_farthest(arbitrary)
    # Second BFS: get diameter path
    far_node2, diameter, parent = bfs_farthest(far_node)
    longest_chain_length = diameter + 1
    
    # Require that the longest chain contains at least 75% of all C atoms.
    if longest_chain_length < 0.75 * n_carbons:
        return False, (f"Carbon chain is too branched: longest chain has {longest_chain_length} atoms "
                       f"but molecule has {n_carbons} carbons")
    
    # Reconstruct the longest chain path.
    path = []
    cur = far_node2
    while cur is not None:
        path.append(cur)
        cur = parent[cur]
    path = list(reversed(path))  # from far_node to far_node2

    # Check that the terminal aldehyde carbon is at one end of the chain.
    if terminal_aldehyde_idx not in (path[0], path[-1]):
        return False, "Terminal aldehyde group is not at the end of the longest carbon chain"
    
    # 9. Evaluate the unsaturation (degree of sp² hybridization) along the carbon chain (excluding the aldehyde carbon).
    # We compute the fraction of sp² carbons in the chain (except the terminal aldehyde carbon).
    chain_without_aldehyde = [idx for idx in path if idx != terminal_aldehyde_idx]
    if chain_without_aldehyde:
        sp2_count = 0
        for idx in chain_without_aldehyde:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() == rdchem.HybridizationType.SP2:
                sp2_count += 1
        sp2_fraction = sp2_count / len(chain_without_aldehyde)
        if sp2_fraction > 0.3:
            return False, (f"Chain has high unsaturation (sp² fraction = {sp2_fraction:.2f}) "
                           "which is uncharacteristic of a typical fatty aldehyde")
    
    return True, "Molecule qualifies as a fatty aldehyde: acyclic, mainly C/H, with a long chain and a single terminal aldehyde group"

# Example test cases (uncomment below to run tests):
# test_smiles = [
#     "[H]C(=CC=O)C(O)CCCCC",    # 4-hydroxynon-2-enal (should be True)
#     "O=CCCCCCCCCC/C=C\\CCCCCCCC",  # 11Z-Eicosenal (should be True)
#     "CCCCCCCCCCCC=O",         # dodecanal (should be True)
#     "O=C/C(=C\\CCCCCC/C=C\\C/C=C\\CCCCC)/CCCCCCCCCCCCC",  # (2Z,10Z,13Z)-2-tridecylnonadeca-2,10,13-trienal (should be True)
# ]
# for s in test_smiles:
#     valid, reason = is_fatty_aldehyde(s)
#     print(f"SMILES: {s}\nResult: {valid}, Reason: {reason}\n")