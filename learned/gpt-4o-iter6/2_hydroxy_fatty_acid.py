"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy group in the alpha- or 2-position relative to the carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Check for alpha/2-hydroxy group
    # General pattern for alpha-hydroxy positioning: longest carbon chain ending in CO group, 2nd carbon having OH
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[CX4,cX3][CX4,cX3](O)[CX3](=O)O")
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "No 2-hydroxy group found at the alpha position"

    return True, "Contains 2-hydroxy group with carboxylic acid"

# Example usage
print(is_2_hydroxy_fatty_acid("CCCCCCC[C@@H](O)C(O)=O"))  # (R)-2-hydroxynonanoic acid