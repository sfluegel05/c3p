"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has an ester functional group with a decanoic acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ester functional group pattern
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester functional group found"

    # Decanoate chain pattern
    decanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if not mol.HasSubstructMatch(decanoate_pattern):
        return False, "No decanoate chain found (10-carbon chain)"

    return True, "Contains ester group with a decanoic acid (10-carbon) chain"