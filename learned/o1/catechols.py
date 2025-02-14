"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: Catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is any compound containing an o-diphenol component,
    which is a ring with two hydroxyl groups on adjacent carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains catechol moiety, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define catechol SMARTS pattern
    # Matches two ring-connected carbons each bearing a hydroxyl group
    catechol_smarts = '[C;R]-[C;R]([OH])-[C;R]([OH])'
    catechol_pattern = Chem.MolFromSmarts(catechol_smarts)

    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains catechol moiety (o-diphenol component)"
    else:
        return False, "Does not contain catechol moiety (o-diphenol component)"