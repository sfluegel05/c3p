"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined by the presence of exactly three non-overlapping carboxylic acid groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxylic acid SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    # Search for carboxylic acid groups in the molecule
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acids = len(carboxylic_acid_matches)

    # We need exactly three separate carboxylic acid groups
    if num_carboxylic_acids == 3:
        # Ensure they are independent and not overlapped contextually
        unique_atoms = set()
        for match in carboxylic_acid_matches:
            unique_atoms.update(match)

        # Unique atoms should be three times the count of groups if they are all distinct
        if len(unique_atoms) == num_carboxylic_acids * 3:
            return True, "Contains exactly three non-overlapping carboxylic acid groups"
        else:
            return False, f"Potential overlap detected in carboxylic acid groups, total distinct atoms: {len(unique_atoms)}"

    return False, f"Contains {num_carboxylic_acids} carboxylic acid groups, expected exactly three"