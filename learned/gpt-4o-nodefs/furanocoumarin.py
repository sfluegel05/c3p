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

    # Refined Furanocoumarin core patterns
    core_patterns = [
        Chem.MolFromSmarts("O=c1cc2oc3ccccc3cc2oc1"),  # Recognizing common benzofuran-2-one structure
        Chem.MolFromSmarts("O=c1oc2ccccc2cc1=O")        # Coumarin-like structure fused to additional ring
    ]

    # Check each pattern
    for pattern in core_patterns:
        if pattern is None:
            continue  # Skip invalid patterns
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a core structure typical of furanocoumarins"
    
    return False, "Does not contain key furanocoumarin structures"