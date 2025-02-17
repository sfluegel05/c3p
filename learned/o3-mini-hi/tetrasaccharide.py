"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: Tetrasaccharide
Definition:
  A tetrasaccharide is an oligosaccharide comprising four monosaccharide units.
  Here we assume that each monosaccharide unit appears as a cyclic structure (typically a furanose or pyranose)
  having one ring oxygen (atomic number 8) and the remaining atoms typically carbons (atomic number 6).
  In addition, the four candidate rings must be connected by glycosidic (or bridging) bonds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import numpy as np

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    It checks that there are exactly four candidate monosaccharide rings (with 5- or 6-membered rings
    containing exactly 1 oxygen) and that those rings are connected via glycosidic bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a tetrasaccharide, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Helper function: check if a ring is a candidate monosaccharide ring.
    def is_monosaccharide_ring(atom_indices):
        # Only consider 5-membered (furanose) or 6-membered (pyranose) rings.
        ring_size = len(atom_indices)
        if ring_size not in [5, 6]:
            return False
        atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
        o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if ring_size == 5:
            return (o_count == 1 and c_count == 4)
        if ring_size == 6:
            return (o_count == 1 and c_count == 5)
        return False

    # Collect all rings that satisfy the monosaccharide criteria.
    candidate_rings = []
    for ring in atom_rings:
        if is_monosaccharide_ring(ring):
            candidate_rings.append(ring)
    
    if len(candidate_rings) != 4:
        return False, f"Found {len(candidate_rings)} monosaccharide-like rings; expected 4 for a tetrasaccharide."

    # Compute an interatomic distance matrix.
    # (Note: RDKit returns a numpy array with distances in terms of bond paths.)
    dm = Chem.GetDistanceMatrix(mol)
    
    # Build a graph (as an adjacency list) between candidate rings.
    # We connect two rings if the minimum bond-distance between an atom in one ring and
    # an atom in the other ring is small (below a threshold). Glycosidic bonds are usually
    # very short (often 1-2 bonds apart) so we choose a threshold of 3.
    threshold = 3
    n = len(candidate_rings)
    # graph: dictionary of node index to list of connected nodes.
    graph = {i: [] for i in range(n)}
    
    for i in range(n):
        for j in range(i+1, n):
            # Compute the minimum distance between any two atoms (one from each ring).
            min_dist = np.inf
            for atom_i in candidate_rings[i]:
                for atom_j in candidate_rings[j]:
                    d = dm[atom_i, atom_j]
                    if d < min_dist:
                        min_dist = d
                    # early break if bond directly connected
                    if min_dist <= 1:
                        break
                if min_dist <= 1:
                    break
            if min_dist <= threshold:
                graph[i].append(j)
                graph[j].append(i)
    
    # Check if the graph of candidate rings is fully connected.
    visited = set()
    def dfs(node):
        visited.add(node)
        for neigh in graph[node]:
            if neigh not in visited:
                dfs(neigh)
    
    dfs(0)
    if len(visited) != n:
        return False, "Monosaccharide-like rings are not connected; glycosidic linkages are missing or not contiguous."
    
    return True, "Contains exactly four connected monosaccharide rings typical for tetrasaccharides."

# Example usage (for testing):
if __name__ == '__main__':
    # Example tetrasaccharide SMILES (from one of the true positive examples)
    sample_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@@H](CO)O[C@H]2O[C@@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O"
    result, reason = is_tetrasaccharide(sample_smiles)
    print(result, reason)