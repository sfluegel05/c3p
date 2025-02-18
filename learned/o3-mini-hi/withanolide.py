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
    (typically four fused rings; the nucleus itself should contain 17 carbons) and
    a side chain that forms a lactone (cyclic ester) ring. The total molecule should have at least 28 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a withanolide, False otherwise.
        str: Explanation / reason for the classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check total carbon count; allow substitutions so require at least 28 C atoms.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 28:
        return False, f"Total carbon count is {total_carbons}, which is too low for a typical withanolide."
    
    # Retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # List of tuples (atom indices) for each ring.
    
    # Filter rings that are characteristic of steroid nucleus: 5- or 6-membered rings.
    rings_filtered = [set(ring) for ring in all_rings if len(ring) in {5, 6}]
    if not rings_filtered:
        return False, "No 5- or 6-membered rings detected; expected a steroid nucleus."
    
    # Group fused rings (rings sharing at least one atom) into clusters.
    # We will build clusters of rings (each cluster is a set of atom indices).
    clusters = []
    for ring in rings_filtered:
        found_cluster = None
        for cluster in clusters:
            if cluster & ring:  # if there is any overlap
                found_cluster = cluster
                break
        if found_cluster:
            found_cluster.update(ring)
        else:
            clusters.append(set(ring))
    
    # Heuristically assume the largest fused cluster is the steroid nucleus.
    nucleus_atoms = max(clusters, key=lambda s: len(s))
    
    # Additionally, we expect the steroid nucleus to arise from about 4 fused rings.
    # Count how many rings in rings_filtered are largely contained (>=80%) in the nucleus.
    fused_ring_count = sum(1 for ring in rings_filtered if ring.issubset(nucleus_atoms))
    if fused_ring_count < 4:
        return False, f"Only {fused_ring_count} fused rings detected in nucleus; expected at least 4 for a steroid nucleus."
    
    # Count how many carbons are in the nucleus (core steroid nucleus is expected to have 17 carbons).
    nucleus_carbons = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if nucleus_carbons != 17:
        return False, f"Steroid nucleus has {nucleus_carbons} carbons, expected 17."

    # Look for a lactone ring, i.e. a cyclic ester.
    # We use a SMARTS that matches an ester where the oxygen is part of a ring (r5 or r6).
    lactone_smarts = "[$([CX3](=O)[OX2;r5]),$([CX3](=O)[OX2;r6])]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    
    lactone_found = False
    # For each lactone candidate, check if the ester fragment is not completely part of the steroid nucleus.
    for match in lactone_matches:
        match_set = set(match)
        if not match_set.issubset(nucleus_atoms):
            lactone_found = True
            break
    
    if not lactone_found:
        return False, "No side chain lactone detected. Withanolides require a modified side chain that forms a lactone ring."
    
    # If all conditions are met, classify the molecule as a withanolide.
    return True, "Molecule is a withanolide: Contains a steroid nucleus with 17 core carbons and a side chain lactone."

# Example usage:
# result, reason = is_withanolide("C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1")
# print(result, reason)