"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin typically has a steroid backbone with hydroxyl groups and attached sugar moieties.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a steroid saponin, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a SMARTS pattern to identify the steroid backbone
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C1)CCC4C3CC(O)CC4C2")
    if steroid_pattern is None:
        return None, None  # Safety check if pattern creation failed
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Use a SMARTS pattern to identify sugar moieties (e.g., glucopyranose ring)
    sugar_pattern = Chem.MolFromSmarts("C1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moieties attached"
    
    # Generalize pattern for additional hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts("O[C;X2]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Hydroxysteroid backbone lacks additional hydroxyl groups"
    
    return True, "Contains a hydroxysteroid backbone with attached sugar moieties consistent with a steroid saponin"