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
    A withanolide is defined as a C28 steroid lactone with a fused steroid nucleus 
    (typically four fused rings) and a lactone formed via a modified side chain.
    
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
    
    # Count carbons: withanolides are defined as C28 steroids.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 28:
        return False, f"Carbon count is {carbon_count}, expected 28 for a typical withanolide."
    
    # Retrieve ring information from the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each entry is a tuple of atom indices forming a ring
    if len(rings) < 4:
        return False, f"Found {len(rings)} rings, but a steroid nucleus should contain at least 4 fused rings."
    
    # Check for a fused ring system: at least one pair of rings must share atoms
    fused_found = False
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            if set(rings[i]).intersection(rings[j]):
                fused_found = True
                break
        if fused_found:
            break
    if not fused_found:
        return False, "No fused ring system detected; expected fused rings indicative of a steroid nucleus."
    
    # Look for a lactone ring.
    # We search for an ester fragment, using a SMARTS pattern for a carbonyl attached to an oxygen.
    # Then confirm that the matching atoms are all part of the same ring (making it cyclic, i.e. lactone).
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    matches = mol.GetSubstructMatches(ester_pattern)
    lactone_found = False
    for match in matches:
        # For each found ester fragment, see if all its atoms belong to any one ring.
        for ring in rings:
            if all(idx in ring for idx in match):
                lactone_found = True
                break
        if lactone_found:
            break
            
    if not lactone_found:
        return False, "No lactone ring detected. Withanolides require a side chain lactone functionality."
    
    # If all conditions are met, class the molecule as a withanolide.
    return True, "Molecule is a withanolide: C28 steroid with a fused steroid nucleus and a lactone ring."
    
# Example usage (the SMILES strings provided can be tested here if desired):
# result, reason = is_withanolide("C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1")
# print(result, reason)