"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol contains exactly two hydroxy groups (-OH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Check for exactly two hydroxyl groups
    if len(hydroxyl_matches) == 2:
        return True, "Contains exactly two hydroxyl groups"
    else:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need exactly 2"