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

    # More flexible steroid core patterns that allow for variations
    steroid_core_patterns = [
        # Basic steroid core with flexible bond types and atom types
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        # Alternative pattern for different steroid variations
        "[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]4~[#6]3~[#6]2~[#6]1"
    ]

    has_steroid_core = False
    for pattern in steroid_core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_steroid_core = True
            break

    if not has_steroid_core:
        return False, "No steroid core structure found"

    # More comprehensive 3β-hydroxy patterns
    beta_hydroxy_patterns = [
        # General 3β-OH pattern with stereochemistry
        '[C@@H;$(C[C@]1[C@@H]2)|$(C[C@@]1[C@H]2)]([OH1])',
        # Alternative pattern for different numbering
        '[C@H;$(C[C@@]1[C@H]2)|$(C[C@]1[C@@H]2)]([OH1])',
        # Pattern allowing for different representations
        '[CH1;R]([OH1])[CH2][CH2][C;R]',
        # Pattern for cyclic systems with beta-OH
        '[C@@;R]([OH1])([#6;R])[#6;R]'
    ]

    has_3beta_hydroxy = False
    for pattern in beta_hydroxy_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_3beta_hydroxy = True
            break

    if not has_3beta_hydroxy:
        return False, "No 3beta-hydroxy group found"

    # Structural validation
    ri = mol.GetRingInfo()
    
    # Check for 4 rings minimum (tetracyclic system)
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    # Check for fused ring system
    fused_rings = False
    for ring in ri.AtomRings():
        for other_ring in ri.AtomRings():
            if ring != other_ring:
                if len(set(ring).intersection(set(other_ring))) >= 2:
                    fused_rings = True
                    break
        if fused_rings:
            break

    if not fused_rings:
        return False, "Missing proper fused ring system"

    # Additional check for proper connectivity
    ring_atoms = set()
    for ring in ri.AtomRings():
        ring_atoms.update(ring)
    if len(ring_atoms) < 16:  # Minimum atoms in steroid core
        return False, "Insufficient ring system"

    return True, "Contains steroid core with 3beta-hydroxy group"