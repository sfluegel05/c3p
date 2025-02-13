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

    # More flexible steroid core pattern that matches both 5-alpha and 5-beta steroids
    # Uses recursive SMARTS to match the four-ring system with more flexibility
    steroid_core = Chem.MolFromSmarts(
        '[C,c]12[C,c][C,c][C,c]3[C,c]([C,c]1)[C,c][C,c]4[C,c][C,c][C,c][C,c]4[C,c][C,c]3[C,c]2'
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Match 3-alpha-hydroxy group in context of the steroid core
    # The pattern looks for the A ring with the 3-alpha-OH group
    alpha_3_hydroxy = Chem.MolFromSmarts(
        '[C,c]12[CH2][CH]([OH])[CH2][C,c][C,c]1[C,c][C,c][C,c]2'
    )
    
    if not mol.HasSubstructMatch(alpha_3_hydroxy, useChirality=True):
        # Try alternative pattern for 5-beta series
        alpha_3_hydroxy_alt = Chem.MolFromSmarts(
            '[C,c]12[CH2][C@H]([OH])[CH2][C,c][C,c]1[C,c][C,c][C,c]2'
        )
        if not mol.HasSubstructMatch(alpha_3_hydroxy_alt, useChirality=True):
            return False, "No 3-alpha hydroxy group found"

    # Validate basic steroid characteristics
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    # Check for characteristic steroid ring sizes
    ring_sizes = sorted([len(r) for r in ring_info.AtomRings()])
    six_membered = sum(1 for size in ring_sizes if size == 6)
    five_membered = sum(1 for size in ring_sizes if size == 5)
    
    if not (six_membered >= 2 and five_membered >= 1):
        return False, "Ring pattern not characteristic of steroids"

    # Additional validation for reasonable molecular weight
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if not (250 < mol_weight < 1000):
        return False, "Molecular weight outside typical steroid range"

    # Check for reasonable number of oxygens (at least one for the 3-OH group)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "No oxygen atoms found"

    return True, "Contains steroid core with 3-alpha hydroxy group"