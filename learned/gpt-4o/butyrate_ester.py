"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is defined as any carboxylic ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define comprehensive butyrate ester pattern (butyric acid as ester)
    butyrate_patterns = [
        Chem.MolFromSmarts("CCCC(=O)O[*]"),  # Butanoate ester group
        Chem.MolFromSmarts("CCC(C)C(=O)O[*]"),  # Variants with branching
    ]

    matches = any(mol.HasSubstructMatch(p) for p in butyrate_patterns)
    if not matches:
        return False, "No butyrate ester group found"

    return True, "Contains a butyrate ester group"