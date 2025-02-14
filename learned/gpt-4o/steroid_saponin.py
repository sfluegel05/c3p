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

    # Define SMARTS pattern for a generic steroid backbone
    steroid_pattern = Chem.MolFromSmarts("C1C[C@@H]2[C@@H](C1)CC[C@H]3[C@@H]2CC[C@H]4C3(C)CC[C@]4(C)C")
    
    # Define SMARTS pattern for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[C][OH]")  # Hydroxyl group attached to carbon

    # Define SMARTS pattern for glycosidic linkage (oxygen linking two carbons)
    glycoside_pattern = Chem.MolFromSmarts("O[C;R][C;R]")  # Oxygen linking ring carbons
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for hydroxyl groups
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found"

    # Check for glycosidic linkage (indicating sugar presence)
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycoside linkages found"
    
    return True, "Contains hydroxysteroid backbone with glycoside linkages"