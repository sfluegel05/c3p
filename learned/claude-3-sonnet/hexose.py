"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:24895 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde group (aldohexose)
    or a ketone group (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core hexose patterns
    hexose_patterns = [
        # Pyranose forms
        "O1[C@H]([CH2]O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O",  # alpha-D-glucopyranose
        "O1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1[CH2]O",  # beta-D-glucopyranose
        "C1([CH2]O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1",    # general pyranose
        
        # Furanose forms
        "O1[C@H]([CH2]O)[C@H](O)[C@H](O)[C@H]1O",  # furanose with CH2OH
        "C1([CH2]O)O[C@H](O)[C@H](O)[C@H]1O",      # general furanose
        
        # Open chain forms
        "[CH]=O[C@H](O)[C@H](O)[C@H](O)[C@H](O)[CH2]O",  # aldehyde form
        "[CH2]O[C@H](O)[C@H](O)[C@H](O)C(=O)[CH2]O",     # ketone form
        
        # Deoxy variants
        "O1[C@H]([CH2]O)[C@H](O)[C@H](O)[C@H](O)[C@H]1[CH3]",  # 6-deoxy
        "O1[C@H]([CH3])[C@H](O)[C@H](O)[C@H](O)[C@H]1O"        # 1-deoxy
    ]

    # Exclusion patterns
    exclusion_patterns = [
        # Phosphates
        "[OH]P(=O)([OH])[OH]",
        "COP(=O)([OH])[OH]",
        # Complex derivatives
        "NS(=O)(=O)[OH]",  # Sulfamates
        "NC(=O)",          # Amides
        "[NH3+]",          # Amino sugars (when charged)
        "P(=O)([OH])([OH])[OH]"  # Phosphoric acid
    ]

    # Check for hexose core structure
    found_core = False
    for pattern in hexose_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_core = True
            break
    
    if not found_core:
        return False, "No hexose core structure found"

    # Check for excluding groups
    for pattern in exclusion_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            return False, "Contains modifications incompatible with simple hexose"

    # Count carbons in main chain
    carbon_chain_pattern = Chem.MolFromSmarts("[C]-[C]-[C]-[C]-[C]-[C]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Missing required six-carbon chain"

    # Count total carbons (allowing for some substitutions)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 6:
        return False, f"Too few carbons for hexose ({carbon_count})"
    if carbon_count > 7:  # Stricter limit on carbons
        # Check if extra carbons are part of simple methyl substitutions
        methyl_pattern = Chem.MolFromSmarts("C[C]")
        methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
        if carbon_count - methyl_count > 6:
            return False, "Too many non-methyl carbons"

    # Check for minimum oxygen content (allowing for deoxy sugars)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 4:  # Allow for deoxy sugars
        return False, f"Too few oxygens for hexose ({oxygen_count})"

    # Verify presence of either aldehyde or ketone group
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    if not (mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(ketone_pattern)):
        # Check for hemiacetal/hemiketal forms
        hemiacetal_pattern = Chem.MolFromSmarts("O[CH]O")
        if not mol.HasSubstructMatch(hemiacetal_pattern):
            return False, "Missing required carbonyl or hemiacetal group"

    return True, "Valid hexose structure with appropriate substitution pattern"