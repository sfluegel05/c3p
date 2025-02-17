"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python3
"""
Classifies: Carotenoid
Definition: Carotenoids are tetraterpenoids (typically having a ~C40 core) derived from psi,psi‐carotene.
They characteristically contain an extended, non‐aromatic conjugated polyene system.
Only atoms C, H, O (and optionally P) are allowed.
This improved version requires:
  – >= 40 carbon atoms in the whole molecule,
  – a molecular weight of at least 400 Da,
  – and an extended, linear conjugated polyene system of at least 10 connected sp² carbons.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict, deque

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Requirements:
      - Only atoms C, H, O, and optionally P are allowed.
      - Must have a sufficiently large carbon skeleton (>=40 carbons) indicative of a tetraterpenoid core.
      - Must have a molecular weight of at least 400 Da.
      - Must feature an extended, non‐aromatic conjugated polyene system of at least 10 sp² carbons.
    
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

    # Allowed atoms are: hydrogen, carbon, oxygen, and phosphorus.
    allowed_atoms = {1, 6, 8, 15}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains disallowed heteroatom: {atom.GetSymbol()}"

    # Count carbons (this is our first filter: we require at least 40 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 40:
        return False, f"Too few carbon atoms ({carbon_count}) to be a carotenoid (expected ~40 in the core)"

    # Check molecular weight (most carotenoids have MW > 400 Da)
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 400:
        return False, f"Molecular weight too low ({mw:.1f} Da) for a carotenoid"

    # Identify sp2-hybridized carbon atoms (candidates for the polyene chain)
    sp2_carbon_indices = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            sp2_carbon_indices.add(atom.GetIdx())

    # Build a graph among sp2 carbons. Two nodes are connected if they are bonded by a conjugated bond.
    graph = defaultdict(list)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() in sp2_carbon_indices and a2.GetIdx() in sp2_carbon_indices:
            if bond.GetIsConjugated():
                graph[a1.GetIdx()].append(a2.GetIdx())
                graph[a2.GetIdx()].append(a1.GetIdx())
    # Ensure that isolated sp2 carbons get an empty list entry.
    for idx in sp2_carbon_indices:
        if idx not in graph:
            graph[idx] = []

    if not graph:
        return False, "No conjugated sp2 carbon system detected"

    # Helper function: perform a breadth-first search (BFS) in a given component.
    def bfs(start, component_nodes):
        visited = {start}
        queue = deque([(start, 1)])  # distance: count starting node as 1
        max_distance = 1
        farthest_node = start
        while queue:
            current, dist = queue.popleft()
            if dist > max_distance:
                max_distance = dist
                farthest_node = current
            for neighbor in graph[current]:
                if neighbor in component_nodes and neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, dist + 1))
        return farthest_node, max_distance

    # Identify connected components (subsets) of the conjugated sp2-carbon graph.
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

    # Among the components, find the longest chain length (the diameter) for a non‐aromatic polyene.
    longest_chain = 0
    for comp in components:
        # Exclude a component if any atom in it is aromatic (we require a non‐aromatic conjugated polyene)
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in comp):
            continue
        # For a cyclic component (if every node has degree 2 and size>=3), we take its size.
        is_cycle = all(len(graph[n]) == 2 for n in comp) and (len(comp) >= 3)
        if is_cycle:
            chain_length = len(comp)
        else:
            # Use two passes of BFS to approximate the diameter of the tree-like structure.
            start = next(iter(comp))
            far_node, _ = bfs(start, comp)
            _, chain_length = bfs(far_node, comp)
        if chain_length > longest_chain:
            longest_chain = chain_length

    # Require that the longest conjugated polyene chain has at least 10 sp2 carbons.
    if longest_chain < 10:
        return False, f"No extended conjugated polyene system found (longest chain length = {longest_chain})"

    # If all checks pass, we classify the molecule as a carotenoid.
    msg = (f"Found {carbon_count} carbons, MW {mw:.1f} Da, "
           f"and a conjugated chain of {longest_chain} sp2 carbons consistent with carotenoid structure")
    return True, msg

# Example usage (this section is only run when executing as script)
if __name__ == "__main__":
    # Replace the following SMILES with any test candidate.
    test_smiles = "CC(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C)C)C)C)C)C)C"
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)