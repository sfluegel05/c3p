"""
Classifies: CHEBI:74716 withanolide
"""
#!/usr/bin/env python
"""
Classifies: Withanolide – “Any steroid lactone that is a C28 steroid with a modified side chain forming a lactone ring and its substituted derivatives.”
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    Withanolides are defined as steroid lactones that have a fused steroid nucleus
    (typically a cyclopentanoperhydrophenanthrene system with around 17 carbons in the nucleus)
    and a side chain that forms a lactone (cyclic ester) ring. The overall molecule should be roughly C28 (with some tolerance)
    and must not be heavily glycosylated or decorated.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule meets withanolide criteria, False otherwise.
        str: Explanation of the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 28:
        return False, f"Total carbon count is {total_carbons}, which is too low for a typical withanolide."
    if total_carbons > 35:
        return False, f"Total carbon count is {total_carbons}, which is too high – may indicate glycosylated cardiac glycoside."
    
    # Retrieve ring information (only consider rings of size 5 or 6, typical for steroid rings).
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # each is a tuple of atom indices
    rings_filtered = [set(ring) for ring in all_rings if len(ring) in {5, 6}]
    if not rings_filtered:
        return False, "No 5- or 6-membered rings detected; expected a steroid nucleus."
    
    # Group overlapping rings into clusters (fused ring systems).
    clusters = []
    for ring in rings_filtered:
        merged = False
        for cluster in clusters:
            if cluster & ring:  # if they share at least one atom
                cluster.update(ring)
                merged = True
                break
        if not merged:
            clusters.append(set(ring))
    
    # Assume that the largest fused-cluster is the steroid nucleus.
    nucleus_atoms = max(clusters, key=lambda s: len(s))
    
    # Count how many of our filtered rings are (mostly) contained in the nucleus.
    fused_ring_count = sum(1 for ring in rings_filtered if ring.issubset(nucleus_atoms))
    if fused_ring_count < 2:
        return False, f"Only {fused_ring_count} fused ring(s) detected in nucleus; expected at least 2 for a steroid nucleus."
    
    # Count the carbons in the nucleus. The steroid nucleus is expected to be around 17 carbons (allowing a tolerance).
    nucleus_carbons = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if not (15 <= nucleus_carbons <= 18):
        return False, f"Steroid nucleus has {nucleus_carbons} carbons; expected roughly 17."
    
    # Look for a lactone ring (cyclic ester). The SMARTS here looks for a carbonyl (CX3=O) bound to an O that is in a ring.
    lactone_smarts = "[CX3](=O)[OX2;r]"  # matches an ester oxygen inside a ring (r indicates ring membership)
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    
    lactone_found = False
    # For each match, check that the lactone fragment is not completely part of the nucleus.
    for match in lactone_matches:
        match_set = set(match)
        if not match_set.issubset(nucleus_atoms):
            lactone_found = True
            break
    if not lactone_found:
        return False, "No side chain lactone detected. Withanolides require a modified side chain that forms a lactone ring."
    
    return True, "Molecule is a withanolide: Contains a steroid nucleus (~17 carbons) and a side chain lactone within a total of ~28 carbons."

# Example usage:
# result, reason = is_withanolide("C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1")
# print(result, reason)