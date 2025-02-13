"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: germacranolide
A sesquiterpene lactone based on germacrane skeleton
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for lactone ring fused to ring system
    # This pattern looks for a lactone (O=C-O) that's part of a ring system
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1O[C@@H]2[C@@H]1*"), # Common trans-fused lactone
        Chem.MolFromSmarts("O=C1O[C@H]2[C@H]1*"),   # Alternative stereochemistry
        Chem.MolFromSmarts("O=C1OC2C1*"),           # Generic fused pattern
        Chem.MolFromSmarts("O=C1OC([CH2])C1"),      # α-methylene-γ-lactone
        Chem.MolFromSmarts("O=C1OCC1"),             # Basic γ-lactone
    ]
    
    has_lactone = False
    for pat in lactone_patterns:
        if pat is not None and mol.HasSubstructMatch(pat):
            has_lactone = True
            break
            
    if not has_lactone:
        return False, "No suitable lactone ring found"

    # Check for 10-membered ring (germacrane skeleton)
    # More specific pattern for germacrane skeleton with typical substitution points
    germacrane_pattern = Chem.MolFromSmarts("C1CCCCCCCCC1")  # Basic 10-membered ring
    ring_info = mol.GetRingInfo()
    has_10_ring = False
    
    for ring in ring_info.AtomRings():
        if len(ring) == 10:
            # Additional check for carbons in the 10-membered ring
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            carbon_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
            if carbon_count >= 8:  # Most atoms should be carbon in germacrane
                has_10_ring = True
                break
                
    if not has_10_ring:
        return False, "No suitable 10-membered ring found"

    # Count carbons (should be around 15 for sesquiterpene)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 25:  # Allowing more flexibility for substituents
        return False, f"Carbon count {c_count} not consistent with sesquiterpene lactone"

    # Check for oxygen count (should have at least 2 for lactone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for germacranolide"

    # Calculate number of rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 2:  # Should have at least 2 rings (10-membered + lactone)
        return False, "Insufficient number of rings"

    # Check for typical germacranolide characteristics
    characteristics = []
    
    # Check for α-methylene-γ-lactone pattern (common but not required)
    methylene_lactone = Chem.MolFromSmarts("C=C1C(=O)OC1")
    if mol.HasSubstructMatch(methylene_lactone):
        characteristics.append("α-methylene-γ-lactone")

    # Check for oxidation patterns
    hydroxy_pattern = Chem.MolFromSmarts("OH")
    if mol.HasSubstructMatch(hydroxy_pattern):
        characteristics.append("hydroxyl groups")

    ester_pattern = Chem.MolFromSmarts("CC(=O)O")
    if mol.HasSubstructMatch(ester_pattern):
        characteristics.append("ester groups")

    # Build reason string
    reason = "Contains germacrane skeleton (10-membered ring) with fused lactone"
    if characteristics:
        reason += f" and {', '.join(characteristics)}"

    return True, reason