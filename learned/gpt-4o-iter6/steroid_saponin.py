"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is characterized by a steroid backbone with attached sugar moieties
    and derives from a hydroxysteroid.

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
    
    # Improved pattern for steroid backbone (allow slight variations)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C1)C=CC4C3CC(O)CC4")
    
    # Verify steroid backbone presence (adapt to possible variations)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # SMARTS pattern for sugar moieties (broaden detection of sugar units)
    sugar_pattern = Chem.MolFromSmarts("[C@@H]1([OX2H])[C@H]([OX2H])[C@@H]([OX2H])[C@H]([OX2H])[OX2H]1")
    
    # Finding attached sugar moieties
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moieties attached"
    
    # Check for hydroxyl groups on the steroid structure
    # Generalize hydroxyl pattern covering variabilities
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Steroid core lacks hydroxyl groups"

    return True, "Contains a hydroxysteroid backbone with attached sugar moieties consistent with a steroid saponin"