"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: Tetrasaccharide
Definition:
  A tetrasaccharide is an oligosaccharide comprising four monosaccharide units.
  In this approach we assume that each monosaccharide unit appears as a cyclic structure 
  (typically a furanose or pyranose) having 5 or 6 atoms with exactly one ring oxygen 
  (atomic number 8) and the remaining atoms being carbon (atomic number 6). In addition, 
  the four candidate rings must be connected by glycosidic (bridging) bonds. Here we define
  a glycosidic linkage as a connection where an atom from one sugar (typically the anomeric 
  carbon) bonds to an exocyclic oxygen (not part of its ring) which in turn bonds to a carbon 
  of a second sugar.
"""
from rdkit import Chem
import numpy as np

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    It first finds candidate monosaccharide rings (5- or 6-membered rings with exactly one oxygen)
    and then checks if exactly four are found. It also verifies that the candidate rings are connected
    via glycosidic bonds (defined here as a linkage where an atom from one ring bonds to an exocyclic 
    oxygen that in turn bonds to an atom in a different ring).
    
    Args:
        smiles (str): SMILES string of the molecule
        
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
    atom_rings = ring_info.AtomRings()  # tuple of tuples of atom indices

    # Helper function: check if a ring is a candidate monosaccharide ring.
    def is_monosaccharide_ring(atom_indices):
        # Only consider 5-membered (furanose) or 6-membered (pyranose) rings.
        ring_size = len(atom_indices)
        if ring_size not in [5, 6]:
            return False
        atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
        # Count atomic numbers: oxygen==8 and carbon==6.
        o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if ring_size == 5:
            return (o_count == 1 and c_count == 4)
        if ring_size == 6:
            return (o_count == 1 and c_count == 5)
        return False

    # Collect candidate rings.
    candidate_rings = []
    for ring in atom_rings:
        if is_monosaccharide_ring(ring):
            candidate_rings.append(set(ring))  # store as set for quick membership checks
    
    if len(candidate_rings) != 4:
        return False, f"Found {len(candidate_rings)} monosaccharide-like rings; expected 4 for a tetrasaccharide."

    # Build connectivity graph among candidate rings.
    # For each pair of rings we look for a glycosidic linkage defined as: an atom (a) in ring i 
    # is connected to an exocyclic oxygen (X) (i.e. not in ring i) and that oxygen is bonded to an atom (b)
    # that is in ring j.
    n = len(candidate_rings)
    graph = {i: [] for i in range(n)}
    
    for i in range(n):
        for j in range(i+1, n):
            connected = False
            # Iterate over atoms in ring i.
            for a_idx in candidate_rings[i]:
                a_atom = mol.GetAtomWithIdx(a_idx)
                # Look at neighbors of candidate atom a.
                for neigh in a_atom.GetNeighbors():
                    neigh_idx = neigh.GetIdx()
                    # Skip if neighbor is also in ring i (we want an exocyclic attachment)
                    if neigh_idx in candidate_rings[i]:
                        continue
                    # Require the linking atom to be oxygen (glycosidic linker is O typically)
                    if neigh.GetAtomicNum() != 8:
                        continue
                    # Now check if this oxygen links to an atom in ring j.
                    for nb in neigh.GetNeighbors():
                        # Do not consider the bond going back to a_atom.
                        if nb.GetIdx() == a_idx:
                            continue
                        if nb.GetIdx() in candidate_rings[j]:
                            connected = True
                            break
                    if connected:
                        break
                if connected:
                    break
            if connected:
                graph[i].append(j)
                graph[j].append(i)
    
    # Use depth-first search to see if these four rings form a single connected cluster.
    visited = set()
    def dfs(node):
        visited.add(node)
        for neigh in graph[node]:
            if neigh not in visited:
                dfs(neigh)
    
    dfs(0)
    if len(visited) != n:
        return False, "Monosaccharide-like rings are not connected via glycosidic bonds; linkages are missing or not contiguous."
    
    return True, "Contains exactly four connected monosaccharide rings with glycosidic linkages typical for tetrasaccharides."

# Example usage (for testing):
if __name__ == '__main__':
    # A true positive example from the provided list.
    sample_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@@H](CO)O[C@H]2O[C@@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O"
    result, reason = is_tetrasaccharide(sample_smiles)
    print(result, reason)