"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehyde
Definition: An aldehyde formally arising from reduction of the carboxylic acid of its corresponding fatty acid,
with a long acyclic chain that is mainly carbon/hydrogen and a terminal aldehyde group.
Improvements:
  • Only allows elements C, H, and O.
  • Requires an acyclic molecular structure.
  • Demands at least 6 carbon atoms.
  • Requires that the non‐H (heavy) atoms are predominantly carbon (fraction ≥ 0.75).
  • Requires exactly one overall carbonyl group ([CX3]=[OX1]) in the molecule.
  • The aldehyde carbon in that carbonyl must be terminal (attached to exactly one carbon).
  • The molecule’s carbon connectivity should be “chain‐like” – the longest chain of connected carbons must account for at least 75% of the total carbons.
"""

from rdkit import Chem
from collections import deque

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule qualifies as a fatty aldehyde.

    Steps:
      1. Parse the SMILES string.
      2. Confirm that only allowed elements (C, H, O) are present.
      3. Check that the molecule is acyclic.
      4. Verify that there are at least 6 carbon atoms.
      5. Confirm that heavy atoms (non-H) are predominantly carbons (fraction ≥ 0.75).
      6. Ensure that the molecule contains exactly one carbonyl group ([CX3]=[OX1]).
      7. Confirm that the carbonyl appears as a terminal aldehyde (the aldehyde carbon is attached to exactly one carbon and is not in a ring).
      8. Compute the longest carbon–chain (using the carbon–atom connectivity) and require that it constitutes at least 75% of all carbon atoms.
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      (bool, str): A tuple where the boolean indicates acceptance and the string explains the decision.
    """

    # 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check allowed elements: only C (6), H (1), and O (8)
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: {atom.GetSymbol()}"

    # 3. Ensure the molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s) which is not typical for a fatty aldehyde"
    
    # 4. Count the number of carbon atoms.
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
        return False, f"Molecule has a low carbon fraction ({carbon_fraction:.2f}); it likely contains additional functionalities"

    # 6. Check that there is exactly one carbonyl group overall.
    # Carbonyl SMARTS: any C (sp3 or sp2) double-bonded to oxygen.
    carbonyl_smarts = "[CX3]=[OX1]"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) != 1:
        return False, f"Expected exactly one carbonyl group; found {len(carbonyl_matches)}"
    
    # 7. Check for a terminal aldehyde group.
    # We use a SMARTS that requires one hydrogen on the carbonyl carbon: [CX3H1](=O)
    aldehyde_smarts = "[CX3H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    terminal_aldehyde_count = 0
    for match in aldehyde_matches:
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        # Exclude if the carbonyl carbon is in a ring (it should be terminal and open)
        if aldehyde_c.IsInRing():
            continue
        # Count carbon neighbors (ignoring the oxygen in the C=O).
        carbon_neighbors = [nbr for nbr in aldehyde_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_aldehyde_count += 1
    
    if terminal_aldehyde_count != 1:
        return False, f"Expected exactly one terminal aldehyde group; found {terminal_aldehyde_count}"
    
    # 8. Assess the carbon-chain linearity.
    # Build a graph (adjacency list) of the indices of carbon atoms.
    carbon_indices = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    carbon_graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in carbon_indices and end in carbon_indices:
            carbon_graph[begin].append(end)
            carbon_graph[end].append(begin)
    
    # Check that the carbon subgraph is connected.
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

    # Get connected component from an arbitrary carbon.
    comp = bfs(next(iter(carbon_indices)))
    if len(comp) != len(carbon_indices):
        return False, "Carbon atoms are not in a single connected chain (fatty chain expected to be contiguous)"
    
    # Compute the (tree) diameter of the carbon graph.
    # In an acyclic graph (tree), the diameter can be found by two BFS passes.
    def bfs_farthest(start):
        seen = {start: 0}
        q = deque([start])
        farthest_node = start
        max_dist = 0
        while q:
            cur = q.popleft()
            for nei in carbon_graph[cur]:
                if nei not in seen:
                    seen[nei] = seen[cur] + 1
                    if seen[nei] > max_dist:
                        max_dist = seen[nei]
                        farthest_node = nei
                    q.append(nei)
        return farthest_node, max_dist

    # First BFS from any node.
    far_node, _ = bfs_farthest(next(iter(carbon_indices)))
    # Second BFS to get the diameter (distance in bonds; add 1 to get number of atoms on the path).
    _, diameter = bfs_farthest(far_node)
    longest_chain_length = diameter + 1

    # Require that the longest chain accounts for at least 75% of all carbons.
    if longest_chain_length < 0.75 * n_carbons:
        return False, (f"Carbon chain is too branched: longest chain has {longest_chain_length} atoms "
                       f"but molecule has {n_carbons} carbons in total")
    
    return True, "Molecule qualifies as a fatty aldehyde: acyclic, mainly C/H, with a long chain and a single terminal aldehyde group"

# Example test cases (uncomment to run):
# test_smiles = [
#     "[H]C(=CC=O)C(O)CCCCC",    # 4-hydroxynon-2-enal (expected True)
#     "O=CCCCCCCCCC/C=C\\CCCCCCCC",  # 11Z-Eicosenal (expected True)
#     "CCCCCCCCCCCC=O",         # dodecanal (expected True)
#     "O=C/C(=C\\CCCCCC/C=C\\C/C=C\\CCCCC)/CCCCCCCCCCCCC",  # (2Z,10Z,13Z)-2-tridecylnonadeca-2,10,13-trienal (expected True)
# ]
# for s in test_smiles:
#     res, reason = is_fatty_aldehyde(s)
#     print(f"SMILES: {s}\nResult: {res}, Reason: {reason}\n")