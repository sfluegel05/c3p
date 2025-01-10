"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    This class has an oxo group located at the 2-position and a carboxylate group as the anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more generalized SMARTS pattern for a 2-oxo group implying flexibility around the structure
    pattern = Chem.MolFromSmarts("C(=O)C([#6])[#6]*C(=O)[O-]")  # allows for intervening atoms between groups

    if mol.HasSubstructMatch(pattern):
        return True, "Contains 2-oxo group and carboxylate anion in the correct context for a 2-oxo monocarboxylic acid anion"

    return False, "Pattern for 2-oxo monocarboxylic acid anion not matched"