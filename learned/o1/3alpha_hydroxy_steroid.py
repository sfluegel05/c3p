"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:17855 3alpha-hydroxy steroid
"""

from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.

    A 3alpha-hydroxy steroid is a steroid with a hydroxy group at position 3 in the alpha configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Define steroid backbone SMARTS pattern
        steroid_smarts = "C1CCC2C(C1)CCC3C2CCC4C3CCCC4"  # Steroid skeleton
        steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
        if steroid_pattern is None:
            return False, "Invalid steroid SMARTS pattern"

        # Check for steroid backbone
        if not mol.HasSubstructMatch(steroid_pattern):
            return False, "No steroid backbone found"

        # Define SMARTS for 3alpha-hydroxy group with alpha stereochemistry
        hydroxy_smarts = "[C@H](O)[C@@H]1CC[C@H]2C1"  # 3alpha-hydroxy group
        hydroxy_pattern = Chem.MolFromSmarts(hydroxy_smarts)
        if hydroxy_pattern is None:
            return False, "Invalid hydroxy SMARTS pattern"

        # Check for the 3alpha-hydroxy group
        if not mol.HasSubstructMatch(hydroxy_pattern):
            return False, "No 3alpha-hydroxy group in alpha configuration found"

        return True, "Molecule is a 3alpha-hydroxy steroid"

    except Exception as e:
        return False, f"Error occurred: {str(e)}"