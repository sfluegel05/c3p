"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    This involves identifying a 2-oxo group at the 2-position in relation to a carboxylate anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define substructure pattern for a 2-oxo group adjacent to a carboxylate anion
    # Pattern matches a carbonyl carbon (indicative of the 2-oxo position) linked to another carbon containing the carboxylate anion
    # Allow flexibility with how far they can be due to possible ring or chain connections
    pattern = Chem.MolFromSmarts("[#6]-[#6](=O)-[#6]-[#6](=O)[O-]")
    if not mol.HasSubstructMatch(pattern):
        return False, "2-oxo and carboxylate group not adjacent as required"

    # If pattern matches the substructure, classify as 2-oxo monocarboxylic acid anion
    return True, "Contains a 2-oxo group at the 2-position adjacent to a carboxylate anion"


# Run tests with various SMILES examples to ensure correctness