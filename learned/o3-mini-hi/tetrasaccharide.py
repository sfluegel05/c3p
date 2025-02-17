"""
Classifies: CHEBI:50126 tetrasaccharide
"""
#!/usr/bin/env python
"""
Classifies: Tetrasaccharide

Definition:
  A tetrasaccharide (by our definition) is an oligosaccharide comprising exactly four monosaccharide units.
  For our purpose, a monosaccharide candidate is defined as a ring that is either five-membered (furanose)
  or six-membered (pyranose) meeting the following criteria:
    • It is non‐aromatic.
    • It contains exactly one ring oxygen atom (atomic num 8) and the remaining atoms are carbons.
       For a 5-membered (furanose) ring, 1 oxygen & 4 carbons; for a 6-membered (pyranose) ring, 1 oxygen & 5 carbons.
    
  Glycosidic linkages are then defined as follows:
    • If a carbon atom in one candidate ring is bonded to an oxygen that is not part of the same ring,
      and that oxygen (which has degree ≥ 2) is bonded to a carbon that belongs to a different candidate ring,
      we record a linkage.
    
  The molecule is considered a tetrasaccharide if exactly four candidate rings are found and they are all connected
  via at least one such glycosidic linkage.
  
In this improved version we first “collect” candidate rings using stricter criteria
and then, separately, build a connectivity graph between rings using the explicit glycosidic linkage rule.
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule qualifies as a tetrasaccharide based on our criteria.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a tetrasaccharide, else False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get rings from RDKit using SSSR (smallest set of smallest rings)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # A candidate monosaccharide ring must be non‐aromatic, 5 or 6 atoms in size,
    # with exactly one oxygen in the ring (and correspondingly 4 carbons for 5-membered,
    # 5 carbons for 6-membered).
    def is_candidate_ring(atom_indices):
        ring_size = len(atom_indices)
        if ring_size not in (5, 6):
            return False
        atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
        # Reject aromatic rings
        if any(atom.GetIsAromatic() for atom in atoms):
            return False
        o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if ring_size == 5 and not (o_count == 1 and c_count == 4):
            return False
        if ring_size == 6 and not (o_count == 1 and c_count == 5):
            return False
        return True
    
    candidate_rings = []
    for ring in atom_rings:
        if is_candidate_ring(ring):
            candidate_rings.append(set(ring))
            
    if len(candidate_rings) != 4:
        return False, f"Found {len(candidate_rings)} monosaccharide-like rings; expected 4 for a tetrasaccharide."
    
    # Define a helper that, given two candidate rings, determines if there is at least one glycosidic linkage.
    # According to our rule, check every carbon in ring A:
    def glycosidic_linkage(ringA, ringB):
        for a_idx in ringA:
            atomA = mol.GetAtomWithIdx(a_idx)
            # Only consider carbon atoms in the ring
            if atomA.GetAtomicNum() != 6:
                continue
            for nbr in atomA.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # We expect the bridging oxygen to not be part of the same ring.
                if nbr.GetAtomicNum() == 8 and nbr_idx not in ringA and nbr.GetDegree() >= 2:
                    # Now check if this oxygen is attached to any carbon belonging to ringB.
                    for oxy_nbr in nbr.GetNeighbors():
                        if oxy_nbr.GetIdx() == a_idx:
                            continue
                        if oxy_nbr.GetAtomicNum() == 6 and oxy_nbr.GetIdx() in ringB:
                            return True
        return False

    # Build connectivity graph (each candidate ring is a node, link if glycosidic linkage is found)
    n = len(candidate_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if glycosidic_linkage(candidate_rings[i], candidate_rings[j]) or glycosidic_linkage(candidate_rings[j], candidate_rings[i]):
                graph[i].add(j)
                graph[j].add(i)
    
    # Check that the 4 rings form a connected graph.
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
    # A sample true positive SMILES (from alpha-L-Fucp-(1->2)-beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-Glcp)
    sample_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@@H](CO)O[C@H]2O[C@@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O"
    classification, explanation = is_tetrasaccharide(sample_smiles)
    print(classification, explanation)