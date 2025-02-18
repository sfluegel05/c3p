"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python3
"""
Classifies: Carotenoid
Definition: Carotenoids are tetraterpenoids (typically featuring a ~C40 core) derived from psi,psi‐carotene.
They characteristically contain an extended, non‐aromatic conjugated polyene system.
Only atoms C, H, O (and sometimes P) are allowed.
Molecules that do not have a long (>=10) linear conjugated polyene system or that contain extra heteroatoms are rejected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict, deque

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Requirements:
      - Only atoms C, H, O, and optionally P are allowed.
      - It should have a significant (typically around 40) number of carbon atoms in the non-glycosylated core.
      - It must have an extended, linear conjugated polyene system (>=10 sp2 carbons) that is non-aromatic.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a carotenoid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Allowed atomic numbers: H(1), C(6), O(8), P(15)
    allowed_atoms = {1, 6, 8, 15}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains disallowed heteroatom: {atom.GetSymbol()}"
    
    # Count the total number of carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Set a minimum carbon threshold for a carotenoid core (~40 required; allow some margin for modifications)
    if carbon_count < 35:
        return False, f"Too few carbon atoms ({carbon_count}) to be a carotenoid (expected ~40 in the core)"
    
    # Build set of indices for sp2-hybridized carbon atoms (candidate for a polyene system).
    sp2_carbon_indices = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            sp2_carbon_indices.add(atom.GetIdx())
    
    # Build a conjugated graph among sp2 carbons.
    # We add an edge if both atoms are sp2 carbons and the bond is conjugated.
    graph = defaultdict(list)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() in sp2_carbon_indices and a2.GetIdx() in sp2_carbon_indices:
            if bond.GetIsConjugated():
                graph[a1.GetIdx()].append(a2.GetIdx())
                graph[a2.GetIdx()].append(a1.GetIdx())
    
    # Make sure to include isolated sp2 carbons as nodes (may be needed to build connected components).
    for idx in sp2_carbon_indices:
        if idx not in graph:
            graph[idx] = []
    
    if not graph:
        return False, "No conjugated sp2 carbon system detected"
    
    # Helper: breadth-first search to calculate the longest distance (chain length) in a connected set.
    def bfs(start, component_nodes):
        visited = {start}
        queue = deque([(start, 1)])  # distance starting at 1 (the starting node counts)
        farthest_node = start
        max_distance = 1
        parent = {start: None}
        while queue:
            current, dist = queue.popleft()
            if dist > max_distance:
                max_distance = dist
                farthest_node = current
            for neighbor in graph[current]:
                if neighbor in component_nodes and neighbor not in visited:
                    visited.add(neighbor)
                    parent[neighbor] = current
                    queue.append((neighbor, dist + 1))
        return farthest_node, max_distance, parent

    # Identify connected components among the conjugated sp2 carbons.
    components = []
    seen = set()
    for node in graph:
        if node not in seen:
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
    
    # For each connected component, estimate the length of the longest chain (diameter).
    # We also verify that the atoms in that component are not aromatic.
    longest_chain = 0
    for comp in components:
        # Reject the component if any atom in it is flagged as aromatic.
        # (Carotenoid conjugated polyene should not be part of an aromatic ring system.)
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in comp):
            continue
        # Determine longest chain in the component.
        # If the component is a pure cycle (each node has degree 2 and component size>=3),
        # we simply use the size (although note that cyclic carotenoids exist, their polyene chain is typically linear).
        is_cycle = True
        for n in comp:
            if len(graph[n]) != 2:
                is_cycle = False
                break
        if is_cycle and len(comp) >= 3:
            chain_length = len(comp)
        else:
            start = next(iter(comp))
            far_node, _, _ = bfs(start, comp)
            _, chain_length, _ = bfs(far_node, comp)
        if chain_length > longest_chain:
            longest_chain = chain_length
    
    # Require the longest conjugated polyene chain to have at least 10 sp2 carbons.
    if longest_chain < 10:
        return False, f"No extended conjugated polyene system found (longest chain length = {longest_chain})"
    
    return True, (f"Found {carbon_count} carbons and a conjugated chain of "
                  f"{longest_chain} sp2 carbons consistent with carotenoids")

# Example usage:
if __name__ == "__main__":
    # Test with a candidate SMILES string (replace with actual test cases as needed)
    test_smiles = "CC(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C)C)C)C)C)C)C"
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)