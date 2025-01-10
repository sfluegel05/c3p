"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # SMARTS patterns
    # Beta-D-glucuronic acid core with correct stereochemistry
    # [OH]-[CH]-[C](=O)-O where the [CH] is part of the ring
    carboxylic_acid = Chem.MolFromSmarts('[OH1][CH1][C](=O)[OH1]')
    
    # Pyranose ring with beta orientation at C1 and correct stereochemistry
    # Note: The exact stereochemistry pattern is complex, we'll check for basic structure
    pyranose = Chem.MolFromSmarts('O1[C@H][C@H][C@H][C@@H][C@@H]1')
    
    # Check for carboxylic acid group
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
        
    # Check for pyranose ring
    if not mol.HasSubstructMatch(pyranose):
        return False, "No pyranose ring with correct stereochemistry found"
        
    # Count oxygen atoms to ensure we have enough for the glucuronic acid structure
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Insufficient oxygen atoms for glucuronic acid structure"
        
    # Check for glycosidic bond - should have at least one O-C-O pattern
    glycosidic = Chem.MolFromSmarts('O[CH1]O')
    if not mol.HasSubstructMatch(glycosidic):
        return False, "No glycosidic bond found"
        
    # Count hydroxyl groups - should have at least 3 (excluding carboxylic acid)
    hydroxyl = Chem.MolFromSmarts('[OH1][CH1]')
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl))
    if hydroxyl_matches < 3:
        return False, "Insufficient hydroxyl groups for glucuronic acid"
        
    # Additional check for ring size
    ring_info = mol.GetRingInfo()
    has_six_membered_ring = False
    for ring_size in ring_info.RingSizes():
        if ring_size == 6:
            has_six_membered_ring = True
            break
    if not has_six_membered_ring:
        return False, "No six-membered ring found"

    return True, "Contains beta-D-glucuronic acid core with glycosidic bond"