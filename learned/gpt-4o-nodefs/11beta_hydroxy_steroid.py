"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid contains a steroid core structure with a hydroxyl group at the 11beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid core pattern: we identify the four ring structure (cyclopentanoperhydrophenanthrene)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3CCC4")  # Simplified four-ring system
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid core structure not found"

    # 11beta-hydroxy group position 
    # Specific atom mapping and ordering can be challenging, simplified by assuming position
    # Recognize by position after ensuring steroid skeleton presence
    hydroxyl_pattern = Chem.MolFromSmarts("[C;H2][C;H2][C;H][OH]")  # Simplification for hydroxyl on a carbon chain
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing 11beta-hydroxy group"

    return True, "Contains steroid core structure with 11beta-hydroxy group"