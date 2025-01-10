"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol should contain exactly two hydroxy groups, allowing for diverse structural contexts.
    
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

    # Pattern to detect hydroxy groups (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count the number of hydroxy groups
    num_hydroxy_groups = len(hydroxy_matches)

    # A diol should have exactly two -OH groups
    if num_hydroxy_groups != 2:
        return False, f"Contains {num_hydroxy_groups} hydroxy groups, need exactly 2"

    return True, "Contains exactly two hydroxy groups forming a suitable diol"