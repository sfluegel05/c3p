"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    A 2-oxo monocarboxylic acid is characterized by:
    - A carbonyl group (C=O) at the alpha (2nd) position relative to a carboxylic acid.
    - The general structure is like RC(=O)C(=O)OH, where R can be a variety of groups.

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

    # Define a SMARTS pattern for a general 2-oxo monocarboxylic acid
    # C(=O)C(=O)O pattern: where the first carbon is connected to any non-hydrogen R group.
    pattern = Chem.MolFromSmarts("C(=O)C(=O)O")

    # Check if the pattern matches
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a 2-oxo monocarboxylic acid substructure"

    return False, "Does not contain a 2-oxo monocarboxylic acid substructure"