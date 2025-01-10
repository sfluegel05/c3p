"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has a decanoic acid moiety (10-carbon chain) attached to an alcohol through an ester linkage.

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

    # Look for decanoate pattern (10-carbon chain followed by an ester group)
    decanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)O")
    if not mol.HasSubstructMatch(decanoate_pattern):
        return False, "No decanoate ester pattern found"
    
    return True, "Contains decanoate ester structure"