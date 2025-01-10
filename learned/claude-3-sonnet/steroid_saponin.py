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
    
    # More flexible steroid core patterns
    steroid_patterns = [
        # Basic steroid ABCD ring system with flexible bond types and substitutions
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6,#8]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1",
        
        # Pattern allowing for ketones and alcohols in ring positions
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6,#8]~[#6]~[#6]~[#6]~[#6]~4~[#6](~[#8,#6])~[#6]~3~[#6]~[#6]~2~[#6]~1",
        
        # Pattern for spirostan and related structures
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1~[#6]~5~[#8]~[#6]~[#6]~[#6]~5",
        
        # More general fused ring pattern
        "[#6]~1~2~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~2",
        
        # Pattern focusing on B/C/D rings with flexible substitution
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6](~[#8,#6,#1])~[#6](~[#8,#6,#1])~[#6]~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1"
    ]
    
    has_steroid_core = False
    for pattern in steroid_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            has_steroid_core = True
            break
            
    if not has_steroid_core:
        return False, "No steroid core structure found"

    # Enhanced sugar patterns
    sugar_patterns = [
        # Pyranose ring with hydroxyls
        "O1[C]([OH0,1])[C]([OH0,1])[C]([OH0,1])[C]([OH0,1])[C]1",
        
        # O-glycosidic bond with flexibility
        "[#6]-[#8]-[#6;R]1[#8][#6]([#8,#6])[#6][#6][#6][#6]1",
        
        # Pattern for deoxy sugars
        "[#6;R]1[#8][#6]([#6])[#6]([#8,#6])[#6]([#8,#6])[#6]1[#6]",
        
        # Pattern for sugar chains
        "[#6;R]1[#8][#6][#6][#6][#6][#6]1[#8]-[#6;R]2[#8][#6][#6][#6][#6][#6]2"
    ]
    
    sugar_count = 0
    for pattern in sugar_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            matches = len(mol.GetSubstructMatches(pat))
            sugar_count += matches
            
    if sugar_count == 0:
        return False, "No sugar moieties found"

    # Check for hydroxyl groups (both free and glycosylated)
    oh_pattern = Chem.MolFromSmarts("[OX2H1,OX2R0-]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    
    # Check for glycosidic oxygens
    glycosidic_o = Chem.MolFromSmarts("[#6;R]-[#8]-[#6]")
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
    if mol_wt < 400:
        return False, "Molecular weight too low for steroid saponin"

    # Ring count (including sugar rings)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 6:  # At least 4 for steroid core + minimum 2 sugar rings
        return False, "Too few rings for steroid saponin"

    # If we have steroid core, sugars, and sufficient oxygenation
    if sugar_count >= 1 and has_steroid_core and total_o_features >= 3 and ring_count >= 6:
        return True, "Contains steroid core with glycosidic linkages and multiple hydroxyl groups"
    
    return False, "Missing key steroid saponin features"