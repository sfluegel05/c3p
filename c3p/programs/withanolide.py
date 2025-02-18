"""
Classifies: CHEBI:74716 withanolide
"""
#!/usr/bin/env python
"""
Classifies: Withanolide – “Any steroid lactone that is a C28 steroid with a modified side chain forming a lactone ring and its substituted derivatives.”
An improved algorithm uses overall carbon count, fused ring (steroid) detection, and a lactone ring that lies outside the core ring system.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    Withanolides are defined as steroid lactones that have roughly 28 carbons overall,
    and contain a fused ring system (steroid nucleus) plus a lactone (cyclic ester) attached as a side chain.
    
    This improved approach:
      1. Checks that the overall carbon count is in a reasonable range (~27 to 33).
      2. Identifies fused rings (only rings of size 5 or 6 are considered) as candidate steroid nucleus.
      3. Requires that a lactone ring (identified by a SMARTS) exists that is not completely part of the candidate nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule meets withanolide criteria, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the total number of carbons in the entire molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Withanolides are C28 steroids (allowing slight substitution tolerance)
    if total_carbons < 27:
        return False, f"Total carbon count is {total_carbons}, which is too low for a typical withanolide."
    if total_carbons > 33:
        return False, f"Total carbon count is {total_carbons}, which is too high – possible glycosylated or unrelated structure."
    
    # Identify rings of size 5 or 6 (the typical sizes in steroid rings).
    ring_info = mol.GetRingInfo()
    rings = [set(r) for r in ring_info.AtomRings() if len(r) in {5, 6}]
    if not rings:
        return False, "No 5- or 6-membered rings detected; expected a fused ring system in a steroid."
    
    # Group overlapping rings (fused rings) into clusters.
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
    
    # Choose the largest cluster (by number of atoms) as candidate steroid nucleus.
    nucleus = max(clusters, key=lambda s: len(s))
    
    # Require that the nucleus includes at least two rings.
    nucleus_ring_count = sum(1 for ring in rings if ring.issubset(nucleus))
    if nucleus_ring_count < 2:
        return False, f"Fused ring system (nucleus) has only {nucleus_ring_count} ring(s); expected at least 2 for a steroid nucleus."
    
    # (Optional) Here one may check that the nucleus contains a reasonable number of carbons.
    nucleus_carbons = sum(1 for idx in nucleus if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    # While classical steroids have ~17 carbons in the nucleus,
    # some withanolides are modified so we use a tolerant range.
    if nucleus_carbons < 10:
        return False, f"Steroid nucleus has {nucleus_carbons} carbons; seems too low for a typical steroid core."
    # We don't set an upper bound since substitution may add extra carbons.
    
    # Look for a lactone ring.
    # This SMARTS matches a lactone fragment: a carbonyl (CX3=O) adjacent to an oxygen that is in a ring.
    lactone_smarts = "[CX3](=O)[OX2;r]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone substructure detected; withanolides require a side-chain lactone."
    
    # Check if at least one lactone match is not entirely contained within the nucleus.
    lactone_found = False
    for match in lactone_matches:
        match_set = set(match)
        # If not all atoms of the lactone match belong to the fused nucleus, consider it as a side-chain lactone.
        if not match_set.issubset(nucleus):
            lactone_found = True
            break
    if not lactone_found:
        return False, "Lactone substructure found but appears to be part of the fused ring system rather than a side chain."
    
    return True, f"Molecule is a withanolide: Total carbons = {total_carbons}, nucleus carbons = {nucleus_carbons}, and a side‐chain lactone detected."


# Example usage:
if __name__ == "__main__":
    # Test one of the approved examples (withaperuvin B)
    test_smiles = "OC12C(C(O)(CC1)C(O)(C3OC(=O)C(=C(C3)C)C)C)(CCC4C2CC(O)C5(O)C4(C)C(=O)C=CC5O)C"
    result, reason = is_withanolide(test_smiles)
    print(result, reason)