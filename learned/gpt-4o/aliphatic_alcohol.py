"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroxyl group attachment to aliphatic (non-aromatic) carbon
    aliphatic_oh_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2]!@[O]")
    oh_matches = mol.GetSubstructMatches(aliphatic_oh_pattern)
    if not oh_matches:
        return False, "No aliphatic hydroxyl group found"

    return True, "Contains aliphatic hydroxyl group(s) on a non-aromatic carbon chain"