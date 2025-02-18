"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Acetate ester (CHEBI:47608)
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is a carboxylic ester where the acid component is acetic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for acetate ester group: O-C(=O)-CH3
    acetate_ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)C")
    # Check if the pattern matches any part of the molecule
    if mol.HasSubstructMatch(acetate_ester_pattern):
        return True, "Contains acetate ester group (O-C(=O)-CH3)"
    else:
        return False, "No acetate ester group found"