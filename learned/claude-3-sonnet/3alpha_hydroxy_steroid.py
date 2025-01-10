"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: 3alpha-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the alpha-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core patterns to match different variants
    steroid_patterns = [
        # Basic steroid core (saturated)
        "[C,c]1[C,c][C,c]2[C,c]([C,c]1)[C,c][C,c]3[C,c]([C,c]2)[C,c][C,c]4[C,c][C,c][C,c][C,c]4[C,c]3",
        # Core with possible double bonds
        "[C,c]1[C,c][C,c]2[C,c]([C,c]1)[C,c][C,c]3[C,c]([C,c]2)[C,c][C,c]4[C,c,C=C][C,c][C,c][C,c]4[C,c]3"
    ]
    
    found_core = False
    for pattern in steroid_patterns:
        core_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(core_pattern):
            found_core = True
            break
    
    if not found_core:
        return False, "No steroid core found"

    # Check for 3-alpha-hydroxy group with multiple possible representations
    alpha_3_hydroxy_patterns = [
        # Standard alpha hydroxyl at position 3
        "[OH][C@H]1[CH2][CH2][C]2",
        # Alternative representation
        "[OH][C@@H]1[CH2][CH2][C]2",
        # More general pattern for alpha hydroxyl
        "[OH][C@H,C@@H]1[CH2][CH2][C]"
    ]
    
    found_alpha_hydroxy = False
    for pattern in alpha_3_hydroxy_patterns:
        hydroxy_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(hydroxy_pattern):
            found_alpha_hydroxy = True
            break
            
    if not found_alpha_hydroxy:
        return False, "No 3-alpha hydroxy group found"

    # Additional verification of ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient ring count for steroid structure"

    # Verify the presence of the characteristic four-ring system
    # Look for rings of appropriate size (6,6,6,5 or similar patterns)
    ring_sizes = [len(r) for r in ring_info.AtomRings()]
    if not any(size in [5,6] for size in ring_sizes):
        return False, "Ring sizes not characteristic of steroid structure"

    return True, "Contains steroid core with 3-alpha hydroxy group"