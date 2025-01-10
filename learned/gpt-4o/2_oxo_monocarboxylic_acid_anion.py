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

    # Look for the 2-oxo group pattern (C(=O)) at the 2-position: [C]-[C](=O)-[other]
    oxo_pattern = Chem.MolFromSmarts("[#6][#6](=O)")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 2-oxo group found at the 2-position"

    # Look for a carboxylate anion pattern (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # If both patterns are found, conclude that the molecule is a 2-oxo monocarboxylic acid anion
    return True, "Contains 2-oxo group and carboxylate anion characteristic of 2-oxo monocarboxylic acid anions"