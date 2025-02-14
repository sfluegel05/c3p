"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    A 2-oxo monocarboxylic acid has:
    - A carbonyl group (C=O) at the alpha (2nd) position.
    - A carboxylic acid functional group (C(=O)O).
    - Any R-group can be attached to the remaining positions, given the 2-oxo position is met.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 2-oxo monocarboxylic acids
    # (The central [C] allows for any non-hydrogen atom attached)
    pattern = Chem.MolFromSmarts("[#6][CX3](=O)[CX3](=O)[OX1H]")

    # Check if the pattern matches
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a 2-oxo monocarboxylic acid substructure"

    return False, "Does not contain a 2-oxo monocarboxylic acid substructure"