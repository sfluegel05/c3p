"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python3
"""
Classifies: Carotenoid
Definition: One of a class of tetraterpenoids (C40 core typically) derived from psi,psi‐carotene,
with an extended conjugated polyene chain. (Retinoids, typically C20, or compounds with non‐canonical heteroatoms are excluded.)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are tetraterpenoids (with a ~C40 core) that feature an extended conjugated polyene system.
    They may be cyclized or modified (e.g. by oxidation or glycosylation) but tend to consist of C, H, O
    (and sometimes P in diphosphate derivatives). Retinoids (typically C20) and molecules with other heteroatoms
    (or a disrupted conjugated system) are excluded.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a carotenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Attempt to parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check allowed elements. For typical carotenoids we expect mostly C, H, O
    # and allow P (as in diphosphate derivatives). If any other elements are present,
    # then we do not classify it as a carotenoid.
    allowed_atoms = {1, 6, 8, 15}  # H, C, O, P (atomic numbers)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains disallowed heteroatom: {atom.GetSymbol()}"
    
    # Count number of carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 35:
        return False, f"Too few carbon atoms ({carbon_count}) to be a carotenoid (expected ~40 in the core)"
    
    # --- Determine the longest contiguous conjugated chain ---
    # We will build a graph of carbon atoms that are sp2 hybridized (typical for a polyene)
    # and where the connecting bond is conjugated.
    conjugated_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            # Only consider sp2-hybridized carbons (which will participate in a polyene)
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                conjugated_atoms.add(atom.GetIdx())
    
    # Build an adjacency list for our graph; two sp2 carbons are connected if there is a bond
    # between them and that bond is conjugated.
    graph = {idx: [] for idx in conjugated_atoms}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() in conjugated_atoms and a2.GetIdx() in conjugated_atoms:
            if bond.GetIsConjugated():
                graph[a1.GetIdx()].append(a2.GetIdx())
                graph[a2.GetIdx()].append(a1.GetIdx())
                
    # Helper function to compute the maximum distance (in number of bonds) 
    # from a starting node using BFS.
    def bfs_max_distance(start, nodes_in_component):
        visited = {start}
        queue = deque([(start, 0)])
        max_dist = 0
        while queue:
            cur, dist = queue.popleft()
            if dist > max_dist:
                max_dist = dist
            for nb in graph.get(cur, []):
                if nb in nodes_in_component and nb not in visited:
                    visited.add(nb)
                    queue.append((nb, dist+1))
        return max_dist
    
    # Identify connected components in the conjugated graph.
    seen = set()
    longest_chain_atoms = 0  # count of atoms in the longest conjugated chain
    for node in graph:
        if node in seen:
            continue
        # Build the connected component using BFS
        comp = set()
        queue = deque([node])
        comp.add(node)
        while queue:
            cur = queue.popleft()
            for nb in graph.get(cur, []):
                if nb not in comp:
                    comp.add(nb)
                    queue.append(nb)
        seen.update(comp)
        # For each node in the component, compute maximum distance.
        # The number of atoms in the longest chain = (max number of bonds) + 1.
        for comp_node in comp:
            distance = bfs_max_distance(comp_node, comp)
            chain_length = distance + 1
            if chain_length > longest_chain_atoms:
                longest_chain_atoms = chain_length
    
    # Heuristic: in most carotenoids the extended conjugated polyene covers at least 8 atoms.
    if longest_chain_atoms < 8:
        return False, "No extended conjugated polyene system found (longest chain length = " + str(longest_chain_atoms) + ")"
    
    return True, f"Found {carbon_count} carbon atoms and a conjugated chain of {longest_chain_atoms} sp2 carbons consistent with carotenoids"

# Example usage:
if __name__ == "__main__":
    # Replace this test SMILES string with any carotenoid candidate.
    test_smiles = "CC(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C)C)C)C)C)C)C"
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)