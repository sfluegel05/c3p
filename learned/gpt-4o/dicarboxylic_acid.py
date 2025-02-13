"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid contains exactly two carboxylic acid function groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carboxylic acid (including both protonated and deprotonated forms)
    carboxylic_acid_pattern_protonated = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_pattern_deprotonated = Chem.MolFromSmarts("C(=O)[O-]")

    # Search for matches in the molecule
    protonated_matches = mol.GetSubstructMatches(carboxylic_acid_pattern_protonated)
    deprotonated_matches = mol.GetSubstructMatches(carboxylic_acid_pattern_deprotonated)

    # Calculate total distinct carboxylic acid groups
    total_matches = len(protonated_matches) + len(deprotonated_matches)

    # Check the number of carboxylic acid groups
    if total_matches == 2:
        return True, "Contains exactly two carboxylic acid groups"
    elif total_matches > 2:
        return False, f"Contains {total_matches} carboxylic acid groups, which is more than 2"
    else:
        return False, f"Contains only {total_matches} carboxylic acid groups, fewer than needed"