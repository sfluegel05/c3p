"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: Tetrasaccharide
Definition:
  A tetrasaccharide (by our definition) is an oligosaccharide comprising exactly four monosaccharide units.
  For our purpose, each monosaccharide candidate is assumed to be a 5-membered (furanose) or 6-membered (pyranose)
  ring that contains exactly one ring oxygen (and the accompanying number of carbons), and which presents at least
  one exocyclic oxygen that bridges to another such candidate ring via a C–O–C bond.
  
  The connectivity rule is: if a carbon belonging to one candidate ring is bonded to an oxygen (outside that ring),
  and that oxygen (with degree ≥2) is in turn bonded to a carbon that belongs to a different candidate ring, then
  a glycosidic linkage is defined between those two rings.
  
  Only if exactly four candidate rings are found and they are all connected via at least one such linkage, we
  consider the molecule a tetrasaccharide.
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule qualifies as a tetrasaccharide by the criteria given in the module docstring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a tetrasaccharide, else False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Check if a ring looks like a monosaccharide: 5- or 6-membered,
    # with exactly one oxygen in the ring (and the corresponding number of carbons),
    # and at least one exocyclic oxygen attached to one of its carbon atoms.
    def is_candidate_ring(atom_indices):
        ring_size = len(atom_indices)
        if ring_size not in (5, 6):
            return False
        atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
        o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        # For furanoses: expect 1 O and 4 C; for pyranoses: 1 O and 5 C.
        if ring_size == 5 and not (o_count == 1 and c_count == 4):
            return False
        if ring_size == 6 and not (o_count == 1 and c_count == 5):
            return False
        
        # Check that at least one ring carbon is attached to an oxygen that is external (not in the ring)
        # and that exocyclic oxygen has degree at least 2 (so it connects elsewhere as well).
        external_oxy_found = False
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms as potential anomeric centers.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in atom_indices:
                    continue  # part of the ring
                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() >= 2:
                    # This oxygen, if attached to any other candidate ring, can be a gateway.
                    external_oxy_found = True
                    break
            if external_oxy_found:
                break
        return external_oxy_found

    # Collect candidate rings. Each candidate is represented as a set of its atom indices.
    candidate_rings = []
    for ring in atom_rings:
        if is_candidate_ring(ring):
            candidate_rings.append(set(ring))
    
    if len(candidate_rings) != 4:
        return False, f"Found {len(candidate_rings)} monosaccharide-like rings; expected 4 for a tetrasaccharide."

    # Build connectivity graph among candidate rings
    # We define a glycosidic linkage as a bond in which a carbon atom belonging to one ring is bonded to an oxygen
    # (which is not part of the ring) and that oxygen in turn is bonded to a carbon in a different candidate ring.
    n = len(candidate_rings)
    graph = {i: set() for i in range(n)}
    # Precompute a mapping from atom index to candidate ring indices (there may be overlap rarely)
    atom_to_ring = {}
    for ri, ring in enumerate(candidate_rings):
        for a in ring:
            atom_to_ring.setdefault(a, set()).add(ri)
    
    # Iterate over each candidate ring and each carbon in it.
    for i, ring_i in enumerate(candidate_rings):
        for a_idx in ring_i:
            atom = mol.GetAtomWithIdx(a_idx)
            if atom.GetAtomicNum() != 6:
                continue  # only consider carbon atoms for initiating a linkage
            # Look at neighbors
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # The oxygen candidate must not be part of the same ring.
                if nbr.GetAtomicNum() == 8 and nbr_idx not in ring_i and nbr.GetDegree() >= 2:
                    # Now check the neighbors of this oxygen (excluding our original carbon)
                    for oxy_nbr in nbr.GetNeighbors():
                        if oxy_nbr.GetIdx() == a_idx:
                            continue
                        # If this neighbor is carbon and it belongs to a different candidate ring
                        if oxy_nbr.GetAtomicNum() == 6:
                            for j in atom_to_ring.get(oxy_nbr.GetIdx(), []):
                                if j != i:
                                    graph[i].add(j)
                                    graph[j].add(i)
    # Now check connectivity: all 4 candidate rings should be connected.
    visited = set()
    def dfs(node):
        visited.add(node)
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor)
    dfs(0)
    if len(visited) != n:
        return False, "Monosaccharide-like rings are not connected by glycosidic bonds; linkages are missing or noncontiguous."
        
    return True, "Contains exactly four connected monosaccharide rings with one or more C–O–C glycosidic linkages typical for tetrasaccharides."

# For testing purposes:
if __name__ == '__main__':
    # Use one of the true positive SMILES from the provided outcomes.
    sample_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@@H](CO)O[C@H]2O[C@@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O"
    classification, explanation = is_tetrasaccharide(sample_smiles)
    print(classification, explanation)