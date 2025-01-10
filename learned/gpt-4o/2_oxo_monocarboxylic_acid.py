"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid has a 2-oxo substituent in the alpha position relative to a monocarboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for 2-oxo moiety in alpha position relative to carboxylic acid group
    # Carbon adjacent to carboxyl carbon ['C](=O) representing 2-oxo-keto part: [*][CX3](=O)C(=O)[O]
    two_oxo_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[CX3](=O)O")

    # Check for 2-oxo group in the correct position relative to carboxylic acid
    if not mol.HasSubstructMatch(two_oxo_pattern):
        return False, "No 2-oxo substituent at correct location relative to monocarboxylic acid found"

    return True, "Contains both 2-oxo group and monocarboxylic acid group"