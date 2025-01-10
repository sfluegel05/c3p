"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester specifically incorporates a 10-carbon alkyl chain linked via an ester group.

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
    
    # Check for various forms of decanoate ester patterns
    decanoate_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCC(=O)O"),
        Chem.MolFromSmarts("C(=O)OCCCCCCCCCC"),
        Chem.MolFromSmarts("[R]O[C](=O)CCCCCCCCCC")
    ]
    
    for pattern in decanoate_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains valid decanoate ester structure"
    
    return False, "No decanoate ester pattern found"