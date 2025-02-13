"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid with specific structural characteristics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define some SMARTS patterns for core flavonoid structure:
    flavonoid_pattern = Chem.MolFromSmarts("c1ccc2c(c1)ccc3c2c(c(o3)c4ccc([OX2H])cc4)")
    methoxy_group_pattern = Chem.MolFromSmarts("O[C;H3]c")
    hydroxyl_group_pattern = Chem.MolFromSmarts("[OX2H]")
    sugar_moiety_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@@H](O)[C@H]1O")
    sulfate_group_pattern = Chem.MolFromSmarts("O=S(=O)(O)c")

    # Check for flavonoid core
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Flavonoid core structure not found"

    # Check for hydroxyl groups
    if not mol.HasSubstructMatch(hydroxyl_group_pattern):
        return False, "No hydroxyl groups found, these are essential for anthoxanthins"

    # Methoxy groups are optional but often present; since missing them isn't disqualifying, we won't check forcefully

    # Check for sugar moiety
    sugar_matches = mol.GetSubstructMatches(sugar_moiety_pattern)
    has_sugar = len(sugar_matches) > 0

    # Check for sulfate groups
    sulfate_matches = mol.GetSubstructMatches(sulfate_group_pattern)
    has_sulfate = len(sulfate_matches) > 0

    # Summarize findings
    if has_sugar:
        return True, "Flavonoid core with sugar moiety present"
    elif has_sulfate:
        return True, "Flavonoid core with sulfate groups present"
    else:
        return True, "Flavonoid core with appropriate solubility-affecting groups"

    return False, "Does not match the expected features of an anthoxanthin"