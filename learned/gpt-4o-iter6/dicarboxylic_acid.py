"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two carboxyl groups.

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

    # Carboxylic acid pattern (O=C[O-] or O=C(O) for free carboxylic acids)
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)O")
    counts_of_carboxyl = len(mol.GetSubstructMatches(carboxylate_pattern))

    # Check that we do not mistakenly include esters; ensure free carboxylic acid groups
    carboxylic_free_pattern = Chem.MolFromSmarts("C(=O)O[H]")
    free_matches = len(mol.GetSubstructMatches(carboxylic_free_pattern))
    
    # Confirm true free carboxylic acids by including both carboxylates and their free forms
    total_acid_count = counts_of_carboxyl + free_matches
    
    # Filter dimers or larger acyl chloride misconceptions
    expected_acid_count = 2

    if total_acid_count == expected_acid_count:
        return True, "Contains exactly 2 carboxyl groups, indicating it is a dicarboxylic acid"

    return False, f"Found {total_acid_count} carboxyl groups, need exactly 2"