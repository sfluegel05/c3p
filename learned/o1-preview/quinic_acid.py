"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or its derivative based on its SMILES string.
    Quinic acid is characterized by a cyclohexane ring with hydroxyl or ester groups at positions
    1, 3, 4, and 5, and a carboxylic acid group attached to the ring at position 1.

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

    # Define the SMARTS pattern for quinic acid core
    # This pattern represents a cyclohexane ring with a carboxylic acid at position 1,
    # and hydroxyl or ester groups at positions 3, 4, and 5
    quinic_acid_smarts = """
    [C@@H]1
    (
        [C@@H]
        (
            [C@@H]
            (
                [C@H]
                (
                    [C@H]
                    (
                        [C@H]1
                        [O]
                    )
                    [OX2H0]
                )
                [OX2H0]
            )
            [OX2H0]
        )
        C(=O)[O,H]
    )
    """

    # Remove whitespace and newlines from SMARTS pattern
    quinic_acid_smarts = ''.join(quinic_acid_smarts.split())

    # Create a mol object from the SMARTS pattern
    pattern = Chem.MolFromSmarts(quinic_acid_smarts)

    if not pattern:
        return False, "Invalid SMARTS pattern for quinic acid"

    # Check for substructure match
    if mol.HasSubstructMatch(pattern):
        return True, "Contains quinic acid core structure"
    else:
        return False, "Does not contain quinic acid core structure"