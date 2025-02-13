"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:35757 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing three carboxy groups (-C(=O)OH).

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

    # Define the carboxylic acid group pattern
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    
    # Find all matches of the carboxylic acid group
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Check if there are exactly 3 carboxylic acid groups
    if len(carboxy_matches) == 3:
        return True, "Contains exactly 3 carboxylic acid groups"
    else:
        return False, f"Found {len(carboxy_matches)} carboxylic acid groups, need exactly 3"