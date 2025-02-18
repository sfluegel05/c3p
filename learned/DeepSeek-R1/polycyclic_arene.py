"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: Polycyclic arene (polycyclic aromatic hydrocarbon)
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon)
    based on its SMILES string. A polycyclic arene has multiple fused aromatic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize to detect aromaticity
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL)
    
    # Get all rings in the molecule
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) < 2:
        return False, "Less than 2 rings"
    
    # Identify aromatic rings (all bonds in ring must be aromatic)
    aromatic_rings = []
    for ring in rings:
        is_aromatic = True
        for i in range(len(ring)):
            a1 = ring[i]
            a2 = ring[(i+1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if not bond or not bond.GetIsAromatic():
                is_aromatic = False
                break
        if is_aromatic:
            aromatic_rings.append(ring)
    
    if len(aromatic_rings) < 2:
        return False, f"Found {len(aromatic_rings)} aromatic rings (needs ≥2)"
    
    # Check for fused aromatic rings (shared bond between at least two aromatic rings)
    fused = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        in_rings = 0
        for ar_ring in aromatic_rings:
            if a1 in ar_ring and a2 in ar_ring:
                in_rings += 1
                if in_rings >= 2:
                    fused = True
                    break
        if fused:
            break
    
    if not fused:
        return False, "No fused aromatic rings found"
    
    return True, "Contains multiple fused aromatic rings"