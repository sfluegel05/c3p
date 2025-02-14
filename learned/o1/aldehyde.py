"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: aldehyde
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is a compound with a carbonyl group bonded to one hydrogen atom and one R group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldehyde functional group SMARTS pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[!O]")  # Carbon with one H, double bonded to O, single bonded to any atom except O

    # Check for aldehyde group
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains an aldehyde functional group"
    else:
        return False, "Does not contain an aldehyde functional group"