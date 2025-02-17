"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins – natural pigments whose core structure is built from 
four pyrrole-like rings united through four bridging (methine) carbons to form 
a fused macrocyclic system.
"""

from rdkit import Chem
from itertools import combinations, permutations

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is defined as a natural pigment whose core structure consists of
    four pyrrole-like units (rings that are aromatic, of size 5 or 6 with exactly one N)
    united via four methine (carbon) bridges. In this heuristic we:
      1. Obtain all rings using GetSymmSSSR.
      2. Select candidate rings that are aromatic, of size 5 or 6, and contain exactly one nitrogen.
      3. Build a connectivity (graph) among candidate rings: two rings are “connected” if they share exactly one atom,
         and that shared atom is a carbon.
      4. Search for any 4 candidates that form a cycle (each candidate has exactly two connections within that set).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a porphyrin, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the smallest set of smallest rings.
    sssr = Chem.GetSymmSSSR(mol)
    if not sssr:
        return False, "No rings found in molecule"
    
    # Step 1: Gather candidate pyrrole-like rings.
    candidate_rings = []  # each entry is a set of atom indices
    for ring in sssr:
        ring_atoms = set(ring)
        ring_size = len(ring_atoms)
        # We consider rings of 5 or 6 members; in porphyrin fused systems they can appear as either.
        if ring_size not in (5, 6):
            continue
        
        # Check that the ring is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
            continue
        
        # Count nitrogen atoms in the ring.
        n_count = sum(1 for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_count != 1:
            continue
        
        candidate_rings.append(ring_atoms)
    
    if len(candidate_rings) < 4:
        return False, f"Found only {len(candidate_rings)} candidate pyrrole-like rings; need at least 4 for a porphyrin core."
    
    # Step 2: Build connectivity graph among candidate rings.
    # Two candidate rings are considered connected if they share exactly one atom and that atom is a carbon.
    graph = {i: set() for i in range(len(candidate_rings))}
    for i, j in combinations(range(len(candidate_rings)), 2):
        inter = candidate_rings[i].intersection(candidate_rings[j])
        if len(inter) == 1:
            common_idx = next(iter(inter))
            if mol.GetAtomWithIdx(common_idx).GetAtomicNum() == 6:
                graph[i].add(j)
                graph[j].add(i)
    
    # Step 3: Look for a cycle of 4 candidate rings.
    # A valid cycle of 4 nodes (rings) should have each node connected to exactly two others
    # (when considering only nodes in the set) and the total number of edges is 4.
    from itertools import combinations
    for nodes in combinations(range(len(candidate_rings)), 4):
        # For this 4-node set, check the induced subgraph.
        valid = True
        total_edges = 0
        for n in nodes:
            neighbors = graph[n].intersection(nodes)
            deg = len(neighbors)
            if deg != 2:
                valid = False
                break
            total_edges += deg
        # In an undirected graph, total_edges (counted per node) should be 8 and unique edges = 4.
        if valid and total_edges == 8:
            return True, "Molecule contains 4 fused pyrrole-like rings connected via 4 methine bridges, consistent with a porphyrin core."
            
    return False, "No valid 4-ring cycle (porphyrin core) found in molecule."

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples (the SMILES here is for a generic porphyrin fragment)
    test_smiles = "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3"
    result, reason = is_porphyrins(test_smiles)
    print("Result:", result)
    print("Reason:", reason)