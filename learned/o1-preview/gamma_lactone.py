"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: gamma-lactone
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a lactone having a five-membered lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the gamma-lactone SMARTS pattern
    # Pattern represents a five-membered ring with an ester linkage (-C(=O)-O-)
    gamma_lactone_pattern = Chem.MolFromSmarts('[O;R][C;R](=O)[C;R][C;R][C;R]')
    if gamma_lactone_pattern is None:
        return None, "Could not parse gamma-lactone SMARTS pattern"

    # Check if molecule contains the gamma-lactone pattern
    if mol.HasSubstructMatch(gamma_lactone_pattern):
        return True, "Contains gamma-lactone ring"
    else:
        return False, "Does not contain gamma-lactone ring"