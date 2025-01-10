"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone typically features a ketone group linked to a secondary alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded SMARTS pattern for secondary alpha-hydroxy ketone
    # This pattern attempts to capture more configurations
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[CX4][OX2H]")

    # Check for the secondary alpha-hydroxy ketone pattern
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains secondary alpha-hydroxy ketone structure"
    
    return False, "Does not contain secondary alpha-hydroxy ketone structure"