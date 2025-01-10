"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester contains the tetradecanoate moiety, which is a 14-carbon saturated fatty acid ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Attempt a broader pattern for tetradecanoate substructure
    # Pattern for C14 chain followed by ester linkage 
    # Allow for some flexibility in surrounding context
    tetradecanoate_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)OC"),
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O"),
        Chem.MolFromSmarts("CCCCCCCCCCC(CCCCCO)C(=O)O"),  # Allow variance
    ]

    for pattern in tetradecanoate_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains tetradecanoate moiety"

    return False, "Does not contain tetradecanoate moiety"