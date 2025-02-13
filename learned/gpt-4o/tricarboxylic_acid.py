"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing three carboxy groups (-C(=O)O
    or its anionic form -C(=O)[O-]).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Define patterns for carboxylic acid group and its anionic form
    carboxy_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1-]")
    
    # Find all carboxylic acid groups (protonated and deprotonated)
    carboxy_matches = mol.GetSubstructMatches(carboxy_acid_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Total carboxylic/carboxylate groups
    total_carboxy_groups = len(carboxy_matches) + len(carboxylate_matches)

    # Check if there are exactly three groups (protonated or deprotonated)
    if total_carboxy_groups == 3:
        return True, "Contains exactly three carboxylic acid or carboxylate groups"
    else:
        return False, f"Found {total_carboxy_groups} carboxylic acid or carboxylate groups, need exactly 3"