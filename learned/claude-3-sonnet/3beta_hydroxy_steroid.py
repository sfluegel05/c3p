"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core patterns - more flexible to catch different variations
    steroid_patterns = [
        # Basic steroid core (cyclopentanoperhydrophenanthrene)
        "[CH2,CH]1[CH2,CH][CH2,CH][C]2([CH2,CH]1)[CH2,CH][CH2,CH][C]1([CH2,CH]2)[CH2,CH][CH2,CH][CH2,CH]2[CH2,CH][CH2,CH][CH2,CH][C]12",
        # Alternative pattern allowing for double bonds
        "[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6][#6]4[#6]3[#6]2[#6]1",
        # Pattern for 5α-steroids
        "[CH2,CH]1[CH2,CH][CH2,CH][C@@H]2[CH2,CH][CH2,CH][C@@H]3[CH2,CH][CH2,CH][CH2,CH]4[CH2,CH][CH2,CH][CH2,CH][C@]4(C)[C@H]3[CH2,CH]2[CH2,CH]1",
        # Pattern for 5β-steroids
        "[CH2,CH]1[CH2,CH][CH2,CH][C@H]2[CH2,CH][CH2,CH][C@@H]3[CH2,CH][CH2,CH][CH2,CH]4[CH2,CH][CH2,CH][CH2,CH][C@]4(C)[C@H]3[CH2,CH]2[CH2,CH]1"
    ]

    has_steroid_core = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_steroid_core = True
            break

    if not has_steroid_core:
        return False, "No steroid core structure found"

    # 3β-hydroxy group patterns with explicit stereochemistry
    beta_hydroxy_patterns = [
        # Explicit 3β-OH pattern with correct stereochemistry
        '[C@@H]1([OH1])[CH2][CH2][C@@]2',  # For 5α-steroids
        '[C@@H]1([OH1])[CH2][CH2][C@]2',   # For 5β-steroids
        # Alternative representation
        '[H][C@@]1([OH1])CC[C@@]2',        # For 5α-steroids
        '[H][C@@]1([OH1])CC[C@]2'          # For 5β-steroids
    ]

    has_3beta_hydroxy = False
    for pattern in beta_hydroxy_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_3beta_hydroxy = True
            break

    if not has_3beta_hydroxy:
        return False, "No 3beta-hydroxy group found"

    # Structural validation
    # Count rings
    ri = mol.GetRingInfo()
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    # Additional check for fused ring system
    fused_ring_pattern = Chem.MolFromSmarts("[R2][R2][R2][R2]")
    if not mol.HasSubstructMatch(fused_ring_pattern):
        return False, "Missing proper fused ring system"

    # Check for proper connectivity of rings
    ring_atoms = set()
    for ring in ri.AtomRings():
        ring_atoms.update(ring)
    if len(ring_atoms) < 16:  # Minimum atoms in steroid core
        return False, "Insufficient ring system"

    return True, "Contains steroid core with 3beta-hydroxy group"