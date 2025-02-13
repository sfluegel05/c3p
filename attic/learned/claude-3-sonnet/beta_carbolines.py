"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies compounds as beta-carbolines based on their SMILES structure
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    Beta-carbolines contain a pyrido[3,4-b]indole core structure and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for beta-carboline core structure
    # Pattern 1: Fully aromatic beta-carboline core
    beta_carboline_pattern1 = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3c2nccc3")
    
    # Pattern 2: Partially hydrogenated beta-carboline core (allowing for variations in saturation)
    beta_carboline_pattern2 = Chem.MolFromSmarts("C1CCc2c(C1)[nH]c3c2nccc3")
    beta_carboline_pattern3 = Chem.MolFromSmarts("[#6]1[#6][#6]c2c([#6]1)[nH]c3c2[#7][#6][#6][#6]3")
    
    # Pattern 4: N-substituted variants (common in the examples)
    beta_carboline_pattern4 = Chem.MolFromSmarts("c1ccc2c(c1)n([#6,#1])c3c2nccc3")
    beta_carboline_pattern5 = Chem.MolFromSmarts("[#6]1[#6][#6]c2c([#6]1)n([#6,#1])c3c2[#7][#6][#6][#6]3")

    # Check for any of the patterns
    if mol.HasSubstructMatch(beta_carboline_pattern1):
        return True, "Contains fully aromatic beta-carboline core structure"
    elif mol.HasSubstructMatch(beta_carboline_pattern2):
        return True, "Contains partially hydrogenated beta-carboline core structure"
    elif mol.HasSubstructMatch(beta_carboline_pattern3):
        return True, "Contains modified beta-carboline core structure"
    elif mol.HasSubstructMatch(beta_carboline_pattern4):
        return True, "Contains N-substituted beta-carboline core structure"
    elif mol.HasSubstructMatch(beta_carboline_pattern5):
        return True, "Contains modified N-substituted beta-carboline core structure"

    # Additional checks for identifying potential false negatives
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring systems for beta-carboline structure"

    return False, "Does not contain beta-carboline core structure"