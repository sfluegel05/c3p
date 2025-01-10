"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is classified as a 2-oxo monocarboxylic acid anion
    based on its SMILES string.
    
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

    # Look for carboxylate group pattern (C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Check for presence of adjacent oxo group (C=O) in the 2-position
    oxo_pattern = Chem.MolFromSmarts("[CX3](=O)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No oxo group found"

    # Check if the molecule contains only one carboxylate group
    # and the overall anionic character is due to the carboxylate group
    num_carboxylate = len(mol.GetSubstructMatches(carboxylate_pattern))
    if num_carboxylate != 1:
        return False, f"Contains {num_carboxylate} carboxylate groups, need exactly 1"

    return True, "Molecule has 2-oxo configuration with a monocarboxylate anion"