"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is characterized by a steroid backbone derived from a hydroxysteroid
    with one or more sugar moieties attached (glycosides).

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

    # Comprehensive steroid nucleus pattern (common four-ring steroid structure)
    steroid_nucleus_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C4CCC(C)(C4)C3)C2C1")
    if not mol.HasSubstructMatch(steroid_nucleus_pattern):
        return False, "No steroid nucleus found"
    
    # Check for hydroxyl group(s) on the steroid backbone
    # Broader pattern to capture variations in hydroxyl group placement
    hydroxyl_pattern = Chem.MolFromSmarts("[#6][#8]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count < 1:
        return False, "No hydroxyl group found on the steroid backbone"

    # Check for sugar moieties (glycosidic bonds), more generic pattern
    # Looks for sugars and their attachment points on steroids
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@H]([C@H](O)[C@H](O1)[C@H]O)CO")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moieties (glycosides) found"

    return True, "Molecule contains steroid nucleus with hydroxysteroid component and glycoside attachments"