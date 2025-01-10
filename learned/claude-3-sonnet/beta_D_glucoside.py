"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule contains a beta-D-glucoside moiety.
    Beta-D-glucoside has specific stereochemistry:
    - Beta configuration at anomeric carbon (C1)
    - D-configuration for the sugar
    - Pyranose form (6-membered ring)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core pattern for beta-D-glucopyranoside
    # Using multiple SMARTS patterns to catch different representations
    glucose_patterns = [
        # Pattern 1: Basic beta-D-glucose core with any substituents
        "[OX2,OH0][C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
        # Pattern 2: More flexible pattern allowing substituted oxygens
        "[OX2,OH0][C@@H]1[C@H]([OH0,OH1])[C@@H]([OH0,OH1])[C@H]([OH0,OH1])[C@@H](CO)O1",
        # Pattern 3: Alternative representation
        "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1"
    ]
    
    found_match = False
    for pattern in glucose_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt, useChirality=True):
            found_match = True
            break
    
    if not found_match:
        return False, "No beta-D-glucose moiety found with correct stereochemistry"

    # Essential validation checks:
    
    # 1. Verify pyranose ring
    pyranose = Chem.MolFromSmarts("O1CCCCC1")
    if not mol.HasSubstructMatch(pyranose):
        return False, "No pyranose ring found"

    # 2. Check for beta configuration at anomeric carbon
    # Multiple patterns to catch different representations
    beta_patterns = [
        "[O][C@@H]1O[C@H]",  # Common pattern
        "[O][C@@H]1O[C@@H]", # Alternative pattern
        "O[C@@H]1O[C@H]"     # Another representation
    ]
    
    beta_found = False
    for pattern in beta_patterns:
        beta_config = Chem.MolFromSmarts(pattern)
        if beta_config is not None and mol.HasSubstructMatch(beta_config, useChirality=True):
            beta_found = True
            break
            
    if not beta_found:
        return False, "Anomeric carbon does not have beta configuration"

    # 3. Basic check for presence of required functional groups
    required_groups = [
        "CO",      # Primary alcohol
        "O1CCCCC1" # Pyranose ring
    ]
    
    for group in required_groups:
        pattern = Chem.MolFromSmarts(group)
        if pattern is None or not mol.HasSubstructMatch(pattern):
            return False, f"Missing required group: {group}"

    return True, "Contains beta-D-glucose moiety with correct stereochemistry and beta glycosidic bond"