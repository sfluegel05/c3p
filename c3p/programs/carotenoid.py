"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python3
"""
Classifies: Carotenoid
Definition: Carotenoids are tetraterpenoids (typically around a C40 core) derived from psi,psi‐carotene,
featuring an extended conjugated polyene system. They may be cyclized or modified (by oxidation, glycosylation, etc.) 
but are composed mostly of C, H, O (and sometimes P). Non‐canonical heteroatoms and molecules with too short 
a conjugated chain (or too few carbons) should NOT be classified as carotenoids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Carotenoids are tetraterpenoids (with a typical ~C40 carbon scaffold) that feature an extended, continuous 
    conjugated polyene chain. They normally only contain C, H, O and (occasionally) P. Molecules that do not meet 
    these criteria (for example retinoids, smaller systems, or compounds with extra heteroatoms) are rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a carotenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Allowed atomic numbers: H (1), C (6), O (8), P (15)
    allowed_atoms = {1, 6, 8, 15}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains disallowed heteroatom: {atom.GetSymbol()}"
    
    # Count number of carbons (must be around 40; we use 35 as minimal cutoff)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 35:
        return False, f"Too few carbon atoms ({carbon_count}) to be a carotenoid (expected ~40 in the core)"
    
    # Build a graph of sp2-hybridized carbon atoms that are linked by conjugated bonds.
    # We only consider carbons that are sp2 so that they likely participate in a polyene.
    sp2_carbon_indices = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            sp2_carbon_indices.add(atom.GetIdx())
    
    # Build a dictionary graph: an edge exists if both atoms are sp2 carbons and the bond between them is conjugated.
    graph = defaultdict(list)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() in sp2_carbon_indices and a2.GetIdx() in sp2_carbon_indices:
            if bond.GetIsConjugated():
                graph[a1.GetIdx()].append(a2.GetIdx())
                graph[a2.GetIdx()].append(a1.GetIdx())
    
    # To better estimate the “linear” extent of the conjugated system,
    # we will compute the longest simple (non-repeating) path in the graph.
    # Since the graph is small, we can perform a recursive DFS.
    def dfs(node, visited):
        max_length = 1  # count current node
        for neighbor in graph[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                path_length = 1 + dfs(neighbor, visited)
                if path_length > max_length:
                    max_length = path_length
                visited.remove(neighbor)
        return max_length

    longest_chain = 0
    # It is possible that the polyene systems are disconnected.
    # We try DFS from all nodes that are in our sp2 carbon graph.
    for node in sp2_carbon_indices:
        if node in graph:  # Only if node has conjugated neighbors
            current_length = dfs(node, {node})
            if current_length > longest_chain:
                longest_chain = current_length

    # Require at least 10 sp2 carbons in a continuous (simple) path
    if longest_chain < 10:
        return False, f"No extended conjugated polyene system found (longest chain length = {longest_chain})"
    
    return True, f"Found {carbon_count} carbon atoms and a conjugated chain of {longest_chain} sp2 carbons consistent with carotenoids"

# Example usage:
if __name__ == "__main__":
    # Test with a candidate carotenoid SMILES string (replace with a candidate)
    test_smiles = "CC(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C)C)C)C)C)C)C"
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)