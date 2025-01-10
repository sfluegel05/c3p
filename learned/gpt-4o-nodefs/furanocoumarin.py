"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    Furanocoumarins are a class of organic chemical compounds that consist of 
    a furan ring fused with a coumarin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded Furanocoumarin core patterns
    core_patterns = [
        Chem.MolFromSmarts("O=C1C=CC2=C(O1)C=CC3=CC=CO3C2"),         # Furo[2,3-h]coumarin
        Chem.MolFromSmarts("O=C1C=CC2=C(O1)C=CC3=C(O2)C=CC=C3"),      # Combined furan and coumarin
        Chem.MolFromSmarts("O=C1C=CC2=C(O1)C=CC3=C(O2)C=CC=CO3"),     # Furo[3,2-g]chromen
        Chem.MolFromSmarts("O=C1OC2=CC=3C=C(O)C(=CC3=C2)C=C1")        # Additional coumarin variations
    ]

    # Check each pattern
    for pattern in core_patterns:
        if pattern is None:
            continue  # Skip invalid patterns
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a core structure typical of furanocoumarins"
    
    return False, "Does not contain key furanocoumarin structures"