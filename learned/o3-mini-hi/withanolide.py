"""
Classifies: CHEBI:74716 withanolide
"""
#!/usr/bin/env python
"""
Classifies: Withanolide – “Any steroid lactone that is a C28 steroid with a modified side chain forming 
a lactone ring and its substituted derivatives.”

Improved algorithm:
  1. Enforces an overall total carbon count within a narrow window.
  2. Attempts to detect a fused steroid nucleus by locating overlapping 5- or 6-membered rings.
     (If only 2 such rings are found, we mark the nucleus as borderline.)
  3. Searches for a lactone substructure whose own ring (of size 5 or 6) is not completely part of the nucleus.
  
N.B. This algorithm is heuristic and not perfect—it may return (None, None) if the logic is too ambiguous.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a given SMILES string corresponds to a withanolide.
    
    Withanolides must:
      1. Have a total carbon count roughly equal to C28 (we allow 27-33 carbons for substituted derivatives).
      2. Possess a fused steroid nucleus (ideally at least 3 overlapping 5- or 6-membered rings,
         but borderline cases with 2 fused rings are flagged as suspicious).
      3. Contain a lactone substructure (a cyclic ester) that is at least partly outside the fused nucleus.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule passes our withanolide criteria, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Count total carbons.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 27:
        return False, f"Total carbon count is {total_carbons}, which is too low for a typical withanolide."
    if total_carbons > 33:
        return False, f"Total carbon count is {total_carbons}, which is too high – possible glycosylated or unrelated structure."
    
    # 2. Identify rings of size 5 or 6 (common in steroids)
    ring_info = mol.GetRingInfo()
    # Get all rings (as tuples of atom indices) that have size 5 or 6.
    candidate_rings = [set(ring) for ring in ring_info.AtomRings() if len(ring) in {5, 6}]
    if not candidate_rings:
        return False, "No 5- or 6-membered rings detected; expected a fused ring system in a steroid."
    
    # Cluster overlapping rings into groups (each group is a candidate steroid nucleus).
    clusters = []
    for ring in candidate_rings:
        found_cluster = False
        for cluster in clusters:
            # If the ring overlaps with any ring in the cluster then merge.
            if cluster & ring:
                cluster.update(ring)
                found_cluster = True
                break
        if not found_cluster:
            clusters.append(set(ring))
    
    # For each cluster, count how many candidate rings are fully contained in it.
    cluster_ring_counts = [sum(1 for ring in candidate_rings if ring.issubset(cluster))
                            for cluster in clusters]
    # Choose the cluster with the maximum number of rings as the nucleus.
    max_idx = cluster_ring_counts.index(max(cluster_ring_counts))
    nucleus = clusters[max_idx]
    nucleus_ring_count = cluster_ring_counts[max_idx]
    nucleus_carbons = sum(1 for idx in nucleus if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    if nucleus_ring_count < 3:
        # Borderline: many withanolides should have at least three fused rings.
        return False, f"Fused ring system (nucleus) has only {nucleus_ring_count} ring(s); expected at least 3 for a steroid nucleus."
    
    # 3. Check for a lactone substructure (cyclic ester).
    # We search using SMARTS: a carbonyl (CX3=O) bonded to an oxygen that is in a ring.
    lactone_smarts = "[CX3](=O)[OX2;r]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone substructure detected; withanolides require a side-chain lactone."
    
    # To reduce false positives, we further require that the lactone belongs to a small ring 
    # of typical size (5 or 6) and that at least one atom of the lactone is outside the nucleus.
    lactone_found = False
    ring_sets = [set(ring) for ring in ring_info.AtomRings()]  # all rings
    for match in lactone_matches:
        match_set = set(match)
        # Check if this match is part of some small ring (5 or 6 atoms) [it may be embedded in a larger ring]
        lactone_ring = None
        for r in ring_sets:
            if match_set.issubset(r) and len(r) in {5,6}:
                lactone_ring = r
                break
        if lactone_ring is None:
            continue  # skip lactone matches not clearly part of a typical small ring
        
        # Now require that not all atoms of this lactone ring are in the nucleus.
        if not lactone_ring.issubset(nucleus):
            lactone_found = True
            break
    if not lactone_found:
        return False, "Lactone substructure found but appears fully embedded in the nucleus rather than in a side-chain."
    
    # All criteria passed.
    return True, (f"Molecule is a withanolide: Total carbons = {total_carbons}, "
                  f"nucleus carbons = {nucleus_carbons} (from {nucleus_ring_count} overlapping rings), "
                  "and a side‐chain lactone detected.")

# Optional: if the classification is too ambiguous, one might instead return (None, None).
# For example, if one wishes to flag borderline cases, add:
#   if borderline_condition:
#       return None, None

# Example usage:
if __name__ == "__main__":
    # Test with one of the approved examples (Withaperuvin B)
    test_smiles = "OC12C(C(O)(CC1)C(O)(C3OC(=O)C(=C(C3)C)C)C)(CCC4C2CC(O)C5(O)C4(C)C(=O)C=CC5O)C"
    result, reason = is_withanolide(test_smiles)
    print(result, reason)