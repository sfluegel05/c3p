"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:78334 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmarts

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxyl group (-OH) on the alpha carbon (second position)
    adjacent to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one carboxylic acid group
    carboxylic_pattern = MolFromSmarts("[CX3](=O)[OX2H1]")  # COOH group
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Find alpha carbon (carbon adjacent to carboxylic acid's carbonyl carbon)
    # and check for hydroxyl group
    alpha_hydroxy_pattern = MolFromSmarts("[CX3](=O)-[CH](-[OH])")  # Alpha carbon with hydroxyl
    if