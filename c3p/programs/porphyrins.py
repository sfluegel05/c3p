"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins – natural pigments whose core structure is built from 
four pyrrole-like rings united through four methine (carbon) bridges to form 
a fused macrocyclic system.
"""

from rdkit import Chem
from itertools import combinations

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is defined as a natural pigment whose core structure consists of
    four pyrrole-like units (aromatic rings of size 5 or 6 with exactly one nitrogen)
    united via four methine (carbon) bridges.
    
    This heuristic:
      1. Uses GetSymmSSSR() to obtain all rings.
      2. Selects candidate rings that are aromatic, of size 5 or 6, and contain exactly one nitrogen.
      3. Constructs a connectivity graph among candidate rings: two rings are connected if they share
         exactly one atom, and that atom is a carbon.
      4. Searches for any 4 candidate rings that form a cycle (each ring connected to exactly two others 
         within the set, forming the fused porphyrin core).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a porphyrin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the smallest set of smallest rings (SSSR)
    sssr = Chem.GetSymmSSSR(mol)
    if not sssr:
        return False, "No rings found in molecule"
    
    # Step 1: Gather candidate pyrrole-like rings.
    candidate_rings = []  # each element is a set of atom indices representing a ring
    for ring in sssr:
        ring_atoms = set(ring)
        ring_size = len(ring_atoms)
        # Consider only rings of size 5 or 6.
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
        shared = candidate_rings[i].intersection(candidate_rings[j])
        if len(shared) == 1:
            common_atom = next(iter(shared))
            if mol.GetAtomWithIdx(common_atom).GetAtomicNum() == 6:
                graph[i].add(j)
                graph[j].add(i)
    
    # Step 3: Look for a cycle of 4 candidate rings.
    # In a valid cycle each of the 4 rings should connect to exactly 2 others within the selected set.
    for nodes in combinations(range(len(candidate_rings)), 4):
        valid_cycle = True
        total_connections = 0
        for node in nodes:
            # Count how many connections this node has among the 4 nodes.
            conn = graph[node].intersection(nodes)
            deg = len(conn)
            if deg != 2:
                valid_cycle = False
                break
            total_connections += deg
        # In an undirected graph, total_connections should be 8 (each edge counts twice).
        if valid_cycle and total_connections == 8:
            return True, "Molecule contains a fused 4-ring porphyrin core with bridging methine carbons."
            
    return False, "No valid 4-ring cycle (porphyrin core) found in molecule."

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES (a generic porphyrin fragment)
    test_smiles = "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3"
    result, reason = is_porphyrins(test_smiles)
    print("Result:", result)
    print("Reason:", reason)