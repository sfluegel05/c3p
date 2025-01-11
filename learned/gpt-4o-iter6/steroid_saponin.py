"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is characterized by a steroid backbone with attached sugar moieties
    and derived from a hydroxysteroid.
    
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
    
    # SMARTS pattern for a generic steroid backbone (simplified)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2(CC1)C3C4CCC(C4)C4C(C3)C(C2)C4")
    
    # Verify steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # SMARTS pattern for sugars, such as glucose, (simplified)
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    
    # Find attached sugar moieties
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar moieties attached"
    
    # Check for hydroxyl groups on the steroid core
    hydroxy_pattern = Chem.MolFromSmarts("CC(O)C")  # Simplified, looking for C-C(O)
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Steroid core lacks hydroxyl group"

    return True, "Contains a hydroxysteroid backbone with attached sugar moieties consistent with a steroid saponin"