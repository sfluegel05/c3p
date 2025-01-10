"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion is formed when the carboxy group of a monocarboxylic acid is deprotonated.

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

    # Look for carboxylate anion group [C](=O)[O-]
    carboxylate_anion_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_anion_pattern)
    num_carboxylate = len(carboxylate_matches)

    if num_carboxylate != 1:
        return False, f"Found {num_carboxylate} carboxylate anion groups, need exactly 1"

    # Look for carboxylic acid group [C](=O)O
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acid = len(carboxylic_acid_matches)

    if num_carboxylic_acid > 0:
        return False, f"Found {num_carboxylic_acid} carboxylic acid groups, should have none"

    # Molecule passes the checks
    return True, "Contains exactly one carboxylate anion group derived from a monocarboxylic acid"