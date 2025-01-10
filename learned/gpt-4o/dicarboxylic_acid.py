"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two carboxyl groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")

    # Find matches for carboxylic acid groups
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)

    # Check the number of carboxylic acid groups
    if len(carboxyl_matches) == 2:
        return True, "Contains exactly two carboxyl groups"
    else:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, need exactly 2"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35683',
        'name': 'dicarboxylic acid',
        'definition': 'Any carboxylic acid containing two carboxy groups.'
    }
}