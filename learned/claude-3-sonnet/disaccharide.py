"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:18667 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two monosaccharides joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count > 3:  # Allow up to 3 rings for special cases
        return False, f"Too many rings ({ring_count}), disaccharides typically have 2"
    if ring_count < 1:
        return False, f"Too few rings ({ring_count}), disaccharides need at least 1"

    # Look for sugar ring patterns with more flexible matching
    pyranose = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1]([CR1][OR1]1)")
    furanose = Chem.MolFromSmarts("[CR1]1[CR1][CR1]([CR1][OR1]1)")
    
    sugar_ring_matches = (len(mol.GetSubstructMatches(pyranose)) + 
                         len(mol.GetSubstructMatches(furanose)))
    
    if sugar_ring_matches < 1:
        return False, "No sugar ring patterns found"

    # Check for glycosidic bond patterns - multiple possible patterns
    glycosidic_patterns = [
        Chem.MolFromSmarts("[CR1][OR2][CR1]"),  # Basic glycosidic bond
        Chem.MolFromSmarts("[CR1]O[CR1;!R]"),    # Acyclic linkage
        Chem.MolFromSmarts("[CR1]O[CR0]"),       # Special cases
    ]
    
    has_glycosidic = False
    for pattern in glycosidic_patterns:
        if mol.HasSubstructMatch(pattern):
            has_glycosidic = True
            break
            
    if not has_glycosidic:
        return False, "No glycosidic bond pattern found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8 or c_count > 24:  # Expanded range for various disaccharides
        return False, f"Carbon count ({c_count}) outside range for disaccharides (8-24)"
    
    if o_count < 6 or o_count > 14:  # Adjusted range
        return False, f"Oxygen count ({o_count}) outside range for disaccharides (6-14)"

    # Look for hydroxyl groups - more flexible pattern
    hydroxyl_patterns = [
        Chem.MolFromSmarts("[CX4][OX2H1]"),  # Standard hydroxyl
        Chem.MolFromSmarts("[CR1][OX2H1]"),   # Ring hydroxyl
    ]
    
    total_hydroxyls = 0
    for pattern in hydroxyl_patterns:
        total_hydroxyls += len(mol.GetSubstructMatches(pattern))
    
    if total_hydroxyls < 3:  # Reduced minimum to account for modifications
        return False, f"Too few hydroxyl groups ({total_hydroxyls}), need at least 3"

    # Check for characteristic sugar carbon patterns - multiple patterns
    sugar_patterns = [
        Chem.MolFromSmarts("[CR1]([OR1])([CR1])[CR1]([OR0])"),
        Chem.MolFromSmarts("[CR1]([OR0])([CR1])[CR1]([OR1])"),
        Chem.MolFromSmarts("[CR1]([OR0])([CH2][OR0])[CR1]"),
    ]
    
    sugar_pattern_matches = 0
    for pattern in sugar_patterns:
        sugar_pattern_matches += len(mol.GetSubstructMatches(pattern))
    
    if sugar_pattern_matches < 2:
        return False, "Insufficient sugar structural patterns"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 240 or mol_wt > 800:  # Slightly adjusted range
        return False, f"Molecular weight ({mol_wt}) outside range for disaccharides (240-800)"

    # Additional check for non-sugar atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if n_count > 2 or s_count > 1:  # Allow limited N for amino sugars
        return False, "Too many non-sugar atoms"

    # Final check for overall connectivity
    if ring_count >= 2 and has_glycosidic and total_hydroxyls >= 3:
        return True, "Contains two sugar units connected by glycosidic bond with appropriate hydroxyl groups"
    
    # If we have one ring and enough evidence of a second sugar unit
    if ring_count == 1 and has_glycosidic and total_hydroxyls >= 4 and sugar_pattern_matches >= 3:
        return True, "Contains connected sugar units with one in ring form"

    return False, "Does not match disaccharide pattern"