"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion has a single carboxylate group (C(=O)[O-]).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylate group pattern
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")

    # Identify matches of the carboxylate pattern
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)

    # Ensure there is exactly one carboxylate group
    if len(carboxylate_matches) == 1:
        return True, "Contains exactly one carboxylate group"

    return False, f"Found {len(carboxylate_matches)} carboxylate groups, which is not one"