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

    # Define a more detailed SMARTS pattern for a steroid backbone
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2C1)CC4=C3CCC4")  # Generic steroid structure (cyclopenta[a]phenanthrene)

    # Example SMARTS pattern for a basic sugar moiety (glucopyranose-like)
    sugar_pattern = Chem.MolFromSmarts("CO[C@H]1[C@H](O)C(O)C(O)C1O")  # Glucose-like pyranose

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for the presence of sugars
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)

    if len(sugar_matches) == 0:
        return False, "No sugar moieties found (glycosidic bonds expected)"
    
    return True, "Contains hydroxysteroid backbone with sugar moieties attached"