"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns
    # Carboxylic acid pattern: -C(=O)O
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    # Hydroxy group at the alpha position: C(O)C(=O)O (alpha hydroxy position)
    alpha_hydroxy_pattern = Chem.MolFromSmarts("C(O)C(=O)O")

    # Check for the presence of carboxylic acid group
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for hydroxy group at the alpha position
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "No hydroxy group at 2-position (alpha carbon) next to carboxylic acid"

    return True, "Contains hydroxy group at 2-position with carboxylic acid group"