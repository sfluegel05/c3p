"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    This involves identifying a 2-oxo group at the 2-position and a carboxylate anion where the
    oxo group is strategically placed adjacent to the carboxylate group.

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

    # Define substructure pattern for a carboxylate anion
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate anion found"
    
    # Define substructure pattern for a 2-oxo group
    # Revised pattern to be less restrictive on the backbone and allow additional substitutions
    two_oxo_pattern = Chem.MolFromSmarts("[#6]C(=O)[#6]-[#6](=O)[O-]") 
    if not mol.HasSubstructMatch(two_oxo_pattern):
        return False, "2-oxo group not present at expected position relative to carboxylate"

    # If both patterns match the substructure, classify as 2-oxo monocarboxylic acid anion
    return True, "Contains a 2-oxo group at the 2-position and a carboxylate anion"

# With the defined test cases, user can verify the accuracy in a separate testing environment