"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: CHEBI:24996 octadecadienoic acid
"""

from rdkit import Chem
from collections import deque

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is any straight-chain C18 polyunsaturated fatty acid
    having two C=C double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a straight-chain fatty acid"

    # Build a carbons-only adjacency list
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 18:
        return False, f"Total number of carbon atoms is {len(carbons)}, expected at least 18"

    carbon_indices = [atom.GetIdx() for atom in carbons]
    adjacency = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            adjacency[idx1].append(idx2)
            adjacency[idx2].append(idx1)

    # Check for cycles in the carbons-only graph
    visited = set()
    parent = {}
    def dfs(u, p):
        visited.add(u)
        for v in adjacency[u]:
            if v == p:
                continue
            if v in visited:
                return True  # Cycle detected
            if dfs(v, u):
                return True
        return False

    # Start DFS from any carbon to check for cycles
    has_cycle = dfs(carbon_indices[0], -1)
    if has_cycle:
        return False, "Carbons form a cyclic structure; not a straight-chain fatty acid"

    # Find the longest carbon chain using BFS
    def bfs(start):
        visited = {start}
        queue = deque([(start, 0, [start])])  # (node, distance, path)
        max_dist = 0
        max_node = start
        max_path = [start]
        while queue:
            current, dist, path = queue.popleft()
            if dist > max_dist:
                max_dist = dist
                max_node = current
                max_path = path
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, dist + 1, path + [neighbor]))
        return max_node, max_dist, max_path

    # Perform BFS twice to find the longest path in the acyclic graph
    _, _, path1 = bfs(carbon_indices[0])
    end_node, _, longest_path = bfs(path1[-1])

    # Check if the longest path has 18 carbons
    if len(longest_path) != 18:
        return False, f"Longest carbon chain length is {len(longest_path)}, expected 18"

    # Check for branching (carbons with more than two carbon neighbors)
    for idx in longest_path:
        atom = mol.GetAtomWithIdx(idx)
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if idx == longest_path[0] or idx == longest_path[-1]:
            # Terminal carbons should have only one carbon neighbor
            if len(carbon_neighbors) != 1:
                return False, "Branching detected at terminal carbon; not a straight-chain fatty acid"
        else:
            # Internal carbons should have exactly two carbon neighbors
            if len(carbon_neighbors) != 2:
                return False, "Branching detected; not a straight-chain fatty acid"

    # Count the number of carbon-carbon double bonds along the chain
    double_bonds = 0
    for i in range(len(longest_path) - 1):
        bond = mol.GetBondBetweenAtoms(longest_path[i], longest_path[i+1])
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bonds += 1

    if double_bonds != 2:
        return False, f"Number of C=C double bonds along the chain is {double_bonds}, expected 2"

    # Check for carboxylic acid or carboxylate group at one end
    def has_carboxyl_group(idx):
        atom = mol.GetAtomWithIdx(idx)
        oxygen_double = False
        oxygen_single = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    oxygen_double = True
                elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    oxygen_single = True
        return oxygen_double and oxygen_single

    if not (has_carboxyl_group(longest_path[0]) or has_carboxyl_group(longest_path[-1])):
        return False, "No carboxylic acid group at the end of the chain"

    return True, "Molecule has a straight-chain C18 fatty acid backbone with two C=C double bonds"