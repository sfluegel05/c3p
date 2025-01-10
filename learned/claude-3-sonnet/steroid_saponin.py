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
    
    # More flexible steroid core patterns including common variations
    steroid_patterns = [
        # Basic steroid core with flexible bonds
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1",
        # Spirostan core common in steroid saponins
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1~[#6]~5~[#8]~[#6]~[#6]~[#6]~5",
        # Cholestane-like core
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]"
    ]
    
    has_steroid_core = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_steroid_core = True
            break
            
    if not has_steroid_core:
        return False, "No steroid core structure found"

    # Sugar patterns including various glycosidic linkages
    sugar_patterns = [
        # Pyranose ring with hydroxyl groups
        "[#6;R6]-1-[#6;R6]-[#6;R6]-[#6;R6]-[#6;R6]-[#8;R6]-1",
        # Glycosidic linkage
        "[#6]-[OX2]-[#6;R](-[OX2H1,OX2R])-[#6;R](-[OX2H1,OX2R])",
        # Alternative sugar pattern
        "[#6;R6]1[#8]C([#6;R6])([#6;R6])[#6;R6][#6;R6][#6;R6]1"
    ]
    
    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatches(Chem.MolFromSmarts(pattern)):
            has_sugar = True
            break
            
    if not has_sugar:
        return False, "No sugar moiety found"

    # Count hydroxyl groups (including those in sugars)
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 3:
        return False, "Insufficient hydroxyl groups"

    # Count oxygens (saponins typically have many oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 8:
        return False, "Too few oxygen atoms for a saponin"

    # Check molecular weight (steroid saponins are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for steroid saponin"

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 5:  # Steroid core (4) + at least one sugar ring
        return False, "Too few rings for steroid saponin"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27:  # Steroid core (~21) + sugar (~6)
        return False, "Too few carbons for steroid saponin"

    return True, "Contains steroid core with glycosidic linkages and multiple hydroxyl groups"