"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.

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

    # Acetate ester pattern: CC(=O)O
    acetate_ester_pattern = Chem.MolFromSmarts("CC(=O)O")
    if not mol.HasSubstructMatch(acetate_ester_pattern):
        return False, "No acetate ester group found"

    return True, "Contains acetate ester group"