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

    # Count carbons and oxygens first as quick filter
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8 or c_count > 24:
        return False, f"Carbon count ({c_count}) outside range for disaccharides (8-24)"
    if o_count < 6 or o_count > 14:
        return False, f"Oxygen count ({o_count}) outside range for disaccharides (6-14)"

    # Check for non-sugar atoms
    non_sugar_atoms = sum(1 for atom in mol.GetAtoms() 
                         if atom.GetAtomicNum() not in [1, 6, 8])
    if non_sugar_atoms > 2:  # Allow up to 2 non-CHO atoms (e.g., N for amino sugars)
        return False, "Too many non-sugar atoms"

    # Define comprehensive patterns for sugar rings
    sugar_ring_patterns = [
        # Pyranose patterns
        "[CR1]1[CR1][CR1][CR1]([CR1][OR1]1)",  # Basic pyranose
        "[CR1]1[CR1][CR1][CR1]([CR1][OR1]1)([OR0])",  # With substituent
        # Furanose patterns
        "[CR1]1[CR1][CR1]([CR1][OR1]1)",  # Basic furanose
        "[CR1]1[CR1][CR1]([CR1][OR1]1)([OR0])",  # With substituent
    ]

    # Define comprehensive glycosidic bond patterns
    glycosidic_patterns = [
        # Common 1->4 linkages
        "[CR1]1[OR1][CR1][CR1][CR1]([CR1]1)O[CR1]",
        # 1->6 linkages
        "[CR1][OR1]C[CR1]",
        # 1->2 and 1->3 linkages
        "[CR1]1[OR1][CR1][CR1]([OR2][CR1])[CR1]1",
        # General patterns
        "[CR1][OR2][CR1;!R]",  # Acyclic linkage
        "[CR1]O[CR1]",         # Basic glycosidic
        # Special cases
        "[CR1]O[CR0]",        # Modified sugars
        "[CR1]1O[CR1][CR1]O1"  # Bridged sugars
    ]

    # Check for sugar rings
    sugar_rings = 0
    for pattern in sugar_ring_patterns:
        smarts = Chem.MolFromSmarts(pattern)
        if smarts:
            sugar_rings += len(mol.GetSubstructMatches(smarts))

    if sugar_rings < 1:
        return False, "No sugar ring patterns found"

    # Check for glycosidic bonds
    has_glycosidic = False
    for pattern in glycosidic_patterns:
        smarts = Chem.MolFromSmarts(pattern)
        if smarts and mol.HasSubstructMatch(smarts):
            has_glycosidic = True
            break

    if not has_glycosidic:
        return False, "No glycosidic bond pattern found"

    # Look for hydroxyl groups with more patterns
    hydroxyl_patterns = [
        "[CX4][OX2H1]",    # Standard hydroxyl
        "[CR1][OX2H1]",    # Ring hydroxyl
        "[CR1][OX2][H]",   # Alternative form
        "[CX4][OX2H]",     # Another variation
    ]
    
    total_hydroxyls = 0
    for pattern in hydroxyl_patterns:
        smarts = Chem.MolFromSmarts(pattern)
        if smarts:
            total_hydroxyls += len(mol.GetSubstructMatches(smarts))

    if total_hydroxyls < 3:
        return False, f"Too few hydroxyl groups ({total_hydroxyls}), need at least 3"

    # Check for characteristic sugar carbon patterns
    sugar_carbon_patterns = [
        "[CR1]([OR1])([CR1])[CR1]([OR0])",
        "[CR1]([OR0])([CR1])[CR1]([OR1])",
        "[CR1]([OR0])([CH2][OR0])[CR1]",
        "[CR1]([OR1])([CR1])([CR1][OR0])",
        "[CR1]1[CR1][CR1]([OR0])[CR1][OR1]1"
    ]
    
    sugar_pattern_matches = 0
    for pattern in sugar_carbon_patterns:
        smarts = Chem.MolFromSmarts(pattern)
        if smarts:
            sugar_pattern_matches += len(mol.GetSubstructMatches(smarts))

    # Final classification
    if sugar_pattern_matches >= 2 and has_glycosidic and total_hydroxyls >= 3:
        return True, "Contains two sugar units connected by glycosidic bond with appropriate hydroxyl groups"
    
    # Check molecular weight as final filter
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 240 or mol_wt > 800:
        return False, f"Molecular weight ({mol_wt}) outside range for disaccharides (240-800)"

    return False, "Does not match disaccharide pattern"