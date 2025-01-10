"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: quinic acid
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or its derivative based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid, characterized by a cyclohexane ring with
    hydroxyl groups and a carboxylic acid attached to the ring. Derivatives may have
    esterified hydroxyl groups or modifications on the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is quinic acid or its derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for quinic acid core (cyclohexane ring with hydroxyl groups and carboxylic acid)
    quinic_acid_smarts = """
        [C;R1]1            # Carbon atom in a ring (cyclohexane)
        ([C;R1][C;R1][C;R1][C;R1][C;R1]1)  # Complete the cyclohexane ring
        (
            [O,N]          # Attached to an oxygen or nitrogen (for hydroxyl or ester)
            [C](=O)O       # Carboxylic acid or ester group
        )
        (*)
    """
    # Clean up the SMARTS pattern
    quinic_acid_smarts = ''.join(quinic_acid_smarts.split())

    pattern = Chem.MolFromSmarts(quinic_acid_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for quinic acid core structure
    if mol.HasSubstructMatch(pattern):
        return True, "Contains quinic acid core structure"
    else:
        return False, "Does not match quinic acid core structure"

__metadata__ = {
    'chemical_class': {
        'name': 'quinic acid',
        'definition': 'A cyclitol carboxylic acid.',
    },
}