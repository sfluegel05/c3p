"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # SMARTS pattern for a steroid backbone, typical of a cyclopenta[a]phenanthrene skeleton
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCCC3C4CCC(O)CC4CCC3C2C1")
    
    # Verify steroid backbone presence
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # SMARTS pattern for sugars, e.g., glucose or other pyranoses
    sugar_pattern = Chem.MolFromSmarts("[OX2H][C@H]1[C@H]([OX2H])[C@H]([OX2H])[C@H]([OX2H])[C@@H]1[OX2H]")
    
    # Finding attached sugar moieties
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moieties attached"
    
    # Check for hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Steroid core lacks hydroxyl groups"

    return True, "Contains a hydroxysteroid backbone with attached sugar moieties consistent with a steroid saponin"