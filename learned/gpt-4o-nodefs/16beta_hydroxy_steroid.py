"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid backbone pattern (cyclopenta[a]phenanthrene system)
    steroid_backbone = Chem.MolFromSmarts('C1CC2CCC3C4CCC(CCC4)C3C2CC1')
    
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"
        
    # Detect the 16beta-hydroxy group
    # Pattern reflects possible stereochemistry at C16 position
    # Exact position of hydroxyl requires neighboring analysis
    pattern_16beta_hydroxy = '[C@@H]1(O)C(CCC2C1CCC3C2CC4CCC(C4)C3)C'
    hydroxyl_pattern = Chem.MolFromSmarts(pattern_16beta_hydroxy)
    
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "16beta-hydroxy group not found or incorrectly positioned"

    return True, "Contains steroid backbone with 16beta-hydroxy group"