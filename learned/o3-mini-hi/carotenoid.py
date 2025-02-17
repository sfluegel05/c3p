"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python3
"""
Classifies: Carotenoid
Definition: Carotenoids are tetraterpenoids (typically around a C40 core) derived from psi,psi‐carotene.
They feature an extended conjugated polyene system. Only C, H, O (and sometimes P) are allowed.
Molecules without an extended (>10 sp2 carbons) conjugated system or with extra heteroatoms are rejected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict, deque

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Carotenoids typically have about 40 carbons and an extended, continuous conjugated polyene chain.
    Only atoms C, H, O (and sometimes P) are allowed. Molecules not meeting these criteria are rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a carotenoid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Allowed atomic numbers (H, C, O, P)
    allowed_atoms = {1, 6, 8, 15}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains disallowed heteroatom: {atom.GetSymbol()}"
    
    # Check if there are at least ~35 carbons to meet the typical C40 core requirement.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 35:
        return False, f"Too few carbon atoms ({carbon_count}) to be a carotenoid (expected ~40 in the core)"
    
    # Build a graph of sp2 carbon atoms connected by conjugated bonds.
    # Each node is the index of an sp2 carbon and an edge exists if the bond is conjugated.
    sp2_carbon_indices = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            sp2_carbon_indices.add(atom.GetIdx())
    
    # Build graph (adjacency list) for conjugated sp2 carbons.
    graph = defaultdict(list)
    for bond in mol.GetBonds():
        # Get the indices of the two atoms in the bond.
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() in sp2_carbon_indices and a2.GetIdx() in sp2_carbon_indices:
            if bond.GetIsConjugated():
                graph[a1.GetIdx()].append(a2.GetIdx())
                graph[a2.GetIdx()].append(a1.GetIdx())
    
    # Also include isolated sp2 carbons (with no conjugated neighbor)
    for idx in sp2_carbon_indices:
        if idx not in graph:
            graph[idx] = []
    
    if not graph:
        return False, "No conjugated sp2 carbon system detected"
    
    # Helper function: breadth-first search from a start node within a given set of nodes.
    # Returns a tuple (farthest_node, distance) where distance is measured in number of nodes along the path.
    def bfs(start, component_nodes):
        visited = {start}
        queue = deque([(start, 1)])  # distance counts the starting node as length 1
        farthest_node = start
        max_distance = 1
        while queue:
            current, dist = queue.popleft()
            if dist > max_distance:
                max_distance = dist
                farthest_node = current
            for neighbor in graph[current]:
                # Only traverse nodes within the current connected component
                if neighbor in component_nodes and neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, dist + 1))
        return farthest_node, max_distance

    # Find connected components in the conjugated graph.
    components = []
    seen = set()
    for node in graph:
        if node not in seen:
            # Use BFS/DFS to get component
            comp_nodes = set()
            stack = [node]
            while stack:
                cur = stack.pop()
                if cur not in comp_nodes:
                    comp_nodes.add(cur)
                    for neighbor in graph[cur]:
                        if neighbor not in comp_nodes:
                            stack.append(neighbor)
            seen |= comp_nodes
            components.append(comp_nodes)
    
    # For each component, estimate the maximum length of a simple conjugated chain.
    # For a pure cycle (all nodes degree 2 and component size>=3) we take the chain length as the size of the component.
    # Otherwise, we compute the "diameter" via two BFS passes.
    longest_chain = 0
    for comp in components:
        # Check if the component is a pure cycle.
        is_cycle = True
        for n in comp:
            # In a pure cycle, each node should have exactly 2 neighbors.
            if len(graph[n]) != 2:
                is_cycle = False
                break
        if is_cycle and len(comp) >= 3:
            chain_length = len(comp)
        else:
            # Compute an approximate diameter:
            # Pick an arbitrary node from the component.
            start = next(iter(comp))
            far_node, _ = bfs(start, comp)
            # Then BFS from that farthest node.
            _, chain_length = bfs(far_node, comp)
        if chain_length > longest_chain:
            longest_chain = chain_length

    # Carotenoids are expected to have an extended polyene: at least 10 sp2 carbons in one continuous system.
    if longest_chain < 10:
        return False, f"No extended conjugated polyene system found (longest chain length = {longest_chain})"
    
    return True, f"Found {carbon_count} carbons and a conjugated chain of {longest_chain} sp2 carbons consistent with carotenoids"

# Example usage:
if __name__ == "__main__":
    # Candidate SMILES string for a carotenoid (replace with test cases as needed)
    test_smiles = "CC(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C)C)C)C)C)C)C"
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)