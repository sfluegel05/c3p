"""
Classifies: CHEBI:74716 withanolide
"""
#!/usr/bin/env python
"""
Classifies: Withanolide – “Any steroid lactone that is a C28 steroid with a modified side chain forming a lactone ring and its substituted derivatives.”
An improved algorithm uses overall carbon count, detection of a fused steroid nucleus (requiring at least 3 overlapping 5- or 6-membered rings),
and the presence of a lactone ring that is not entirely embedded in that nucleus.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    
    Withanolides must:
      1. Have an overall carbon count in the expected range (27 to 33 carbons).
      2. Possess a fused ring system (steroid nucleus) composed of at least three overlapping 5- or 6-membered rings.
      3. Contain a lactone ring (cyclic ester) identified by the SMARTS "[CX3](=O)[OX2;r]" that is not completely part of the fused nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule passes withanolide criteria, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Count the total number of carbons.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 27:
        return False, f"Total carbon count is {total_carbons}, which is too low for a typical withanolide."
    if total_carbons > 33:
        return False, f"Total carbon count is {total_carbons}, which is too high – possible glycosylated or unrelated structure."
    
    # 2. Identify rings of sizes 5 or 6 (common in steroids).
    ring_info = mol.GetRingInfo()
    rings = [set(ring) for ring in ring_info.AtomRings() if len(ring) in {5, 6}]
    if not rings:
        return False, "No 5- or 6-membered rings detected; expected a fused ring system in a steroid."
    
    # Group overlapping rings into clusters. Each cluster is a candidate fused nucleus.
    clusters = []
    for ring in rings:
        merged = False
        for cluster in clusters:
            if cluster & ring:
                cluster.update(ring)
                merged = True
                break
        if not merged:
            clusters.append(set(ring))
    
    # For each cluster, count how many of the candidate rings (from our list) lie completely in that cluster.
    cluster_ring_counts = []
    for cluster in clusters:
        count = sum(1 for ring in rings if ring.issubset(cluster))
        cluster_ring_counts.append(count)
    
    # Choose the cluster with the maximum number of fused rings as candidate steroid nucleus.
    max_index = cluster_ring_counts.index(max(cluster_ring_counts))
    nucleus = clusters[max_index]
    nucleus_ring_count = cluster_ring_counts[max_index]
    
    # Require that the nucleus contains at least 3 fused rings.
    if nucleus_ring_count < 3:
        return False, f"Fused ring system (nucleus) has only {nucleus_ring_count} ring(s); expected at least 3 for a steroid nucleus."
    
    # Optionally count the number of carbons in the nucleus.
    nucleus_carbons = sum(1 for idx in nucleus if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # 3. Look for a lactone (cyclic ester) that is not fully contained in the nucleus.
    # This SMARTS finds a carbonyl (CX3=O) bonded to an oxygen that is in a ring.
    lactone_smarts = "[CX3](=O)[OX2;r]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone substructure detected; withanolides require a side-chain lactone."
    
    # At least one lactone match must have a part outside the steroid nucleus.
    lactone_found = False
    for match in lactone_matches:
        match_set = set(match)
        if not match_set.issubset(nucleus):
            lactone_found = True
            break
    if not lactone_found:
        return False, "Lactone substructure found but appears to be part of the fused ring system rather than a side-chain."
    
    return True, (f"Molecule is a withanolide: Total carbons = {total_carbons}, " 
                  f"nucleus carbons = {nucleus_carbons} (from {nucleus_ring_count} rings), "
                  "and a side‐chain lactone detected.")

# Example usage:
if __name__ == "__main__":
    # Test one of the approved examples (withaperuvin B)
    test_smiles = "OC12C(C(O)(CC1)C(O)(C3OC(=O)C(=C(C3)C)C)C)(CCC4C2CC(O)C5(O)C4(C)C(=O)C=CC5O)C"
    result, reason = is_withanolide(test_smiles)
    print(result, reason)