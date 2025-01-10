"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is defined as a compound that contains exactly two hydroxy groups (-OH).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a diol, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Hydroxy groups pattern, limiting to alcohols (non-aromatic/primary)
    hydroxy_pattern = Chem.MolFromSmarts("[C,c;!$(C=O)]O")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count the number of unique hydroxy groups
    num_hydroxy_groups = len(hydroxy_matches)
    
    # A diol should have exactly two hydroxy groups
    if num_hydroxy_groups == 2:
        return True, "Contains exactly two hydroxy groups, classified as diol"
    else:
        return False, f"Found {num_hydroxy_groups} hydroxy groups, diol requires exactly 2"

# Testing the function on provided SMILES strings would follow...