"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is characterized by a hydroxysteroid backbone with attached sugar moieties.

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

    # Define a more flexible SMARTS pattern for a steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[#6]1([#6][#6])[#6][#6]2[#6]3[#6][#6]([#6]1)[#6][#6]4[#6]2[#6]3[#6]([#6]4)[#8]")  # Four-ring steroid-like structure with generic carbon and at least one hydroxyl group

    # SMARTS pattern for sugars (common pyranose ring)
    sugar_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O")  # Simplified pyranose pattern

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for the presence of sugars
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)

    if len(sugar_matches) == 0:
        return False, "No sugar moieties found (glycosidic bonds expected)"
    
    return True, "Contains hydroxysteroid backbone with sugar moieties attached"