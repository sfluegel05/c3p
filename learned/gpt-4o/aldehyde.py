"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a compound contains an aldehyde group (RC(=O)H) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains an aldehyde group, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for an aldehyde group: carbonyl carbon bonded to a hydrogen
    # Improved with an additional check on adjacency to ensure proper context
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#1]")

    # Check if the molecule contains the aldehyde functional group
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains an aldehyde group (RC(=O)H)"
    
    return False, "No aldehyde group found"

# Test cases could be added here for validation