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

    # Define a basic steroid nucleus pattern (perhydrocyclopentanophenanthrene)
    steroid_nucleus_pattern = Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(steroid_nucleus_pattern):
        return False, "No steroid nucleus found"
    
    # Check for hydroxyl group(s) on the steroid backbone
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4;R]O")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found on the steroid backbone"

    # Check for sugar moieties (glycosidic bonds) attached to the structure
    sugar_pattern = Chem.MolFromSmarts("C1OC(CO)C(CO)O1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moieties (glycosides) found"

    # Additional checks can be added if more specific features of steroid saponins are necessary

    return True, "Molecule contains steroid nucleus with hydroxysteroid component and glycoside attachments"