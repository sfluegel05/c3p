"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: Tetrasaccharide
Definition:
  A tetrasaccharide is an oligosaccharide comprising four monosaccharide units.
  In this approach we assume that each monosaccharide unit appears as a cyclic structure 
  (typically a furanose or pyranose) having 5 or 6 atoms with exactly one ring oxygen 
  and that it is connected to at least one other unit by a glycosidic linkage.
  
  A glycosidic linkage is defined here as a bond in which an atom of one candidate ring 
  (usually a carbon) is connected to an exocyclic oxygen (not in that ring) and that oxygen 
  in turn is bonded to an atom that belongs to a different candidate ring.
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule qualifies as a tetrasaccharide.
    It identifies candidate monosaccharide rings (5- or 6-membered rings with exactly one ring oxygen)
    and then filters them further to include only rings that show at least one exocyclic oxygen that could
    participate in a glycosidic bond. It then builds a connectivity graph among those rings using the rule 
    that a candidate glycosidic linkage is present when an atom from ring A is bonded to an oxygen (not part of A)
    that in turn is bonded to an atom in ring B.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a tetrasaccharide, else False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Helper function: checks if a ring (given by atom indices) is a candidate monosaccharide ring.
    # Criteria: ring size = 5 (furanose) or 6 (pyranose), exactly one oxygen in ring,
    # and at least one exocyclic oxygen neighbor attaches to a ring atom.
    def is_candidate_ring(atom_indices):
        ring_size = len(atom_indices)
        if ring_size not in (5, 6):
            return False
        atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
        o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        # For furanoses: 1 O and 4 C; for pyranoses: 1 O and 5 C.
        if ring_size == 5 and not (o_count == 1 and c_count == 4):
            return False
        if ring_size == 6 and not (o_count == 1 and c_count == 5):
            return False
        
        # Check for at least one exocyclic oxygen attached to any ring atom.
        # This oxygen should not belong to the ring.
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in atom_indices:
                    continue
                if nbr.GetAtomicNum() == 8:  # exocyclic oxygen candidate
                    # Optionally, check that the oxygen has at least one neighbor in a ring (will do that later).
                    return True
        return False

    # Collect candidate rings (store as sets of indices)
    candidate_rings = []
    for ring in atom_rings:
        if is_candidate_ring(ring):
            candidate_rings.append(set(ring))
    
    if len(candidate_rings) != 4:
        return False, f"Found {len(candidate_rings)} monosaccharide-like rings; expected 4 for a tetrasaccharide."

    # Build connectivity graph among candidate rings.
    # For each pair of candidate rings, scan atoms in ring_i.
    # For an atom in ring_i, if there is a neighbor that is oxygen (and not in ring_i),
    # and if that oxygen is also bonded to any atom that belongs to ring_j, define a glycosidic linkage.
    n = len(candidate_rings)
    graph = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            found_link = False
            for a_idx in candidate_rings[i]:
                a_atom = mol.GetAtomWithIdx(a_idx)
                # Look at neighbors outside the ring.
                for nbr in a_atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in candidate_rings[i]:
                        continue
                    if nbr.GetAtomicNum() != 8:  # only consider oxygen as linker
                        continue
                    # For this exocyclic oxygen, check if it links to ring j.
                    for o_nbr in nbr.GetNeighbors():
                        # Avoid going back to a_atom.
                        if o_nbr.GetIdx() == a_idx:
                            continue
                        if o_nbr.GetIdx() in candidate_rings[j]:
                            found_link = True
                            break
                    if found_link:
                        break
                if found_link:
                    break
            if found_link:
                graph[i].append(j)
                graph[j].append(i)
    
    # To be a true tetrasaccharide, all 4 candidate rings must be connected via at least one glycosidic bond.
    # Use DFS on the graph.
    visited = set()
    def dfs(node):
        visited.add(node)
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor)
    dfs(0)
    if len(visited) != n:
        return False, "Monosaccharide-like rings are not connected via glycosidic bonds; linkages are missing or not contiguous."
    
    return True, "Contains exactly four connected monosaccharide rings with glycosidic linkages typical for tetrasaccharides."

# For testing purposes:
if __name__ == '__main__':
    # Use one of the true positive SMILES from the provided list.
    sample_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@@H](CO)O[C@H]2O[C@@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O"
    result, reason = is_tetrasaccharide(sample_smiles)
    print(result, reason)