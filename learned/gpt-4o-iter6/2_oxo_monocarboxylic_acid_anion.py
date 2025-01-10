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
    # Pattern matches a carbon doubly bonded to oxygen (oxo group) adjacent to a carbon with a carboxylate group
    pattern = Chem.MolFromSmarts("[#6]-[#6;D3](=O)-[#6](=O)[O-]")
    if not mol.HasSubstructMatch(pattern):
        return False, "2-oxo and carboxylate group not adjacent as expected"

    # If pattern matches the substructure, classify as 2-oxo monocarboxylic acid anion
    return True, "Contains a 2-oxo group at the 2-position adjacent to a carboxylate anion"

# Test this function with SMILES examples to ensure proper classification