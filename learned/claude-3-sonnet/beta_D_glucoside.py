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
    # [OH0] means any O that's not a hydroxyl (glycosidic linkage)
    # Specific stereochemistry at each carbon matching beta-D-glucose
    glucose_core = """
        [OH0][C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](CO)O1
    """
    pattern = Chem.MolFromSmarts(glucose_core.strip())
    if pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # More flexible pattern allowing for substitutions
    glucose_flex = """
        [OH0][C@@H]1[C@H]([OH0,OH1])[C@@H]([OH0,OH1])[C@H]([OH0,OH1])[C@H](CO)O1
    """
    pattern_flex = Chem.MolFromSmarts(glucose_flex.strip())
    if pattern_flex is None:
        return False, "Invalid SMARTS pattern"

    # Check for matches with either pattern
    has_core = mol.HasSubstructMatch(pattern, useChirality=True)
    has_flex = mol.HasSubstructMatch(pattern_flex, useChirality=True)
    
    if not (has_core or has_flex):
        return False, "No beta-D-glucose moiety found with correct stereochemistry"

    # Additional validation checks:
    
    # 1. Verify pyranose ring
    pyranose = Chem.MolFromSmarts("O1CCCCC1")
    if pyranose is None or not mol.HasSubstructMatch(pyranose):
        return False, "No pyranose ring found"

    # 2. Check for CH2OH group
    ch2oh = Chem.MolFromSmarts("CO[H,C,O]")
    if ch2oh is None or not mol.HasSubstructMatch(ch2oh):
        return False, "Missing characteristic CH2OH group"

    # 3. Count oxygens attached to the glucose core
    oxygen_count = Chem.MolFromSmarts("[OX2][CH1,CH2]")
    if oxygen_count is None:
        return False, "Invalid oxygen pattern"
    
    matches = mol.GetSubstructMatches(oxygen_count)
    if len(matches) < 5:  # Need at least 5 oxygens (4 OH + 1 ring O)
        return False, "Insufficient oxygen-containing groups for beta-D-glucose"

    # 4. Verify beta configuration at anomeric carbon
    beta_config = Chem.MolFromSmarts("[O][C@@H]1O[C@H]")
    if beta_config is None or not mol.HasSubstructMatch(beta_config, useChirality=True):
        return False, "Anomeric carbon does not have beta configuration"

    return True, "Contains beta-D-glucose moiety with correct stereochemistry and beta glycosidic bond"