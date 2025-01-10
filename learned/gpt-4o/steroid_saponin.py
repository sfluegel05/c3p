"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Classify a chemical entity as a steroid saponin based on its SMILES string.
    A steroid saponin is defined as any saponin derived from a hydroxysteroid with
    sugar moieties (glycosides).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a steroid saponin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated and more general steroid backbone pattern
    steroid_patterns = [
        Chem.MolFromSmarts("C1CC2CCC3C4CCC(C)(C4)CC3C2C1"),  # Common steroid backbone
        Chem.MolFromSmarts("C1CCC2C(CC3C(C2C1)C4CCC5C(C4)C3)C5")  # Another variant
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No generic steroid nucleus found"

    # Hydroxyl group check on the steroid structure
    hydroxyl_pattern = Chem.MolFromSmarts("C[OH]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count < 1:
        return False, "No hydroxyl groups found on the steroid backbone"

    # Check for any sugar moieties (glycosidic linkages), using broad patterns
    sugar_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@H](CO)[C@H](O)[C@H](O)[C@@H]1O"),  # Basic sugar pattern
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]1")  # Cyclic sugars
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns):
        return False, "No sugar moieties (glycosides) found"
    
    return True, "Molecule is a steroid saponin with steroid nucleus, hydroxysteroid component, and glycoside attachments"