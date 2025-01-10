"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: steroid saponin
Definition: Any saponin derived from a hydroxysteroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic steroid ring patterns - breaking down into smaller substructures
    steroid_patterns = [
        # Basic ABCD ring system with flexibility in bond types
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6][#6]4[#6][#6]3[#6][#6]2[#6]1",
        # Alternative pattern allowing for unsaturation
        "[#6]1[#6][#6]2[#6]=,:[#6][#6]3[#6]=,:[#6][#6]4[#6][#6][#6][#6][#6]4[#6][#6]3[#6][#6]2[#6]1",
        # Pattern for 5α-spirostan type
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6][#6]4[#6][#6]3[#6][#6]2[#6]1[#6]5[#8][#6][#6][#6]5",
        # Pattern focusing on B/C/D rings
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1"
    ]
    
    has_steroid_core = False
    for pattern in steroid_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            has_steroid_core = True
            break
            
    if not has_steroid_core:
        return False, "No steroid core structure found"

    # Sugar and glycosidic linkage patterns
    sugar_patterns = [
        # Pyranose ring
        "O1[C][C][C][C][C]1",
        # O-glycosidic bond
        "[#6]-[#8]-[#6;R]1[#8][#6][#6][#6][#6][#6]1",
        # Alternative sugar pattern with hydroxyls
        "[#6;R]1[#8][#6]([#6][#6][#6][#6]1)([#8])",
        # Specific glycosidic linkage to steroid
        "[#6;R][#8][#6;R]1[#8][#6][#6][#6][#6][#6]1"
    ]
    
    sugar_count = 0
    for pattern in sugar_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            matches = len(mol.GetSubstructMatches(pat))
            sugar_count += matches
            
    if sugar_count == 0:
        return False, "No sugar moieties found"

    # Check for hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    
    # Check for glycosidic oxygens
    glycosidic_o = Chem.MolFromSmarts("[#6;R]-[#8]-[#6;!R]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_o)) if glycosidic_o else 0
    
    total_o_features = oh_matches + glycosidic_matches
    if total_o_features < 3:
        return False, "Insufficient hydroxyl/glycosidic groups"

    # Count total oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for a saponin"

    # Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Lowered threshold to catch smaller steroid saponins
        return False, "Molecular weight too low for steroid saponin"

    # Ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 5:
        return False, "Too few rings for steroid saponin"

    # Additional check for characteristic steroid saponin features
    if sugar_count >= 1 and has_steroid_core and total_o_features >= 3:
        return True, "Contains steroid core with glycosidic linkages and multiple hydroxyl groups"
    
    return False, "Missing key steroid saponin features"