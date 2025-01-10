"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined by the presence of exactly three distinct and non-overlapping 
    carboxylic acid groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxylic acid SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")

    # Search for carboxylic acid groups in the molecule
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    # Initialize set for unique carboxylic acids
    unique_atom_sets = []

    for match in carboxylic_acid_matches:
        unique_atom_set = set(match)

        # Check for overlap with existing groups
        if not any(unique_atom_set.intersection(other_set) for other_set in unique_atom_sets):
            unique_atom_sets.append(unique_atom_set)

    num_distinct_carboxylic_acids = len(unique_atom_sets)

    if num_distinct_carboxylic_acids == 3:
        return True, "Contains exactly three distinct and non-overlapping carboxylic acid groups"
    else:
        return False, f"Found {num_distinct_carboxylic_acids} distinct carboxylic acid groups, expected exactly three"