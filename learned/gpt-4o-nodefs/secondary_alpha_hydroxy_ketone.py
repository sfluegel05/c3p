"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone is defined by a ketone group with an adjacent secondary alcohol group.

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

    # SMARTS pattern for secondary alpha-hydroxy ketone
    # [C](-[OH])(-[C])C(=O)C where [C] has exactly two other carbon neighbors (secondary)
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[C](O)([C])[C](=O)[C]")
    
    # Check for the secondary alpha-hydroxy ketone pattern
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains secondary alpha-hydroxy ketone structure"

    return False, "Does not contain secondary alpha-hydroxy ketone structure"