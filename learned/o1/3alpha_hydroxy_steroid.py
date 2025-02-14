"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:17855 3alpha-hydroxy steroid
"""

from rdkit import Chem
from rdkit.Chem import rdqueries

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

        # Define a general steroid backbone SMARTS pattern
        # Steroids have a cyclopentanoperhydrophenanthrene skeleton: three fused six-membered rings and one five-membered ring fused together
        steroid_smarts = "[$([R2]1[R2][R2][R2][R2][R2][R2][R2]1)]"  # A fused ring system
        steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCCC4")  # Typical steroid skeleton
        if not mol.HasSubstructMatch(steroid_pattern):
            # Try a more general pattern
            steroid_pattern = Chem.MolFromSmarts("[$([R2]1[R2][R2][R2][R2][R2][R2][R2]1)]")
            if not mol.HasSubstructMatch(steroid_pattern):
                return False, "No steroid backbone found"

        # Define SMARTS for 3alpha-hydroxy group with alpha stereochemistry
        # Carbon at position 3 with hydroxyl group in alpha configuration
        hydroxy_smarts = "[C@@H](O)[C@]1([H])CC[C@H]2C1"
        hydroxy_pattern = Chem.MolFromSmarts(hydroxy_smarts)
        if hydroxy_pattern is None:
            return False, "Invalid hydroxy SMARTS pattern"

        # Check for the 3alpha-hydroxy group
        if not mol.HasSubstructMatch(hydroxy_pattern):
            return False, "No 3alpha-hydroxy group in alpha configuration found"

        return True, "Molecule is a 3alpha-hydroxy steroid"

    except Exception as e:
        return False, f"Error occurred: {str(e)}"