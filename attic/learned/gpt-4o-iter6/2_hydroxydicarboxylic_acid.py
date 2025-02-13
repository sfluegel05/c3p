"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid must have two carboxylic acid groups and a hydroxy group
    on the alpha carbon to one of the carboxylic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the two carboxylic acid groups
    carboxy_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    if len(carboxy_matches) < 2:
        return False, "Less than two carboxylic acid groups found"

    # Look for the hydroxy group on the alpha carbon
    hydroxyl_alpha_carboxy_pattern = Chem.MolFromSmarts("C(O)C(=O)O")
    alpha_oh_matches = mol.GetSubstructMatches(hydroxyl_alpha_carboxy_pattern)
    if not alpha_oh_matches:
        return False, "No hydroxy group found on the alpha carbon to a carboxylic acid group"
    
    return True, "Contains two carboxylic acid groups and a hydroxy group on the alpha carbon"