"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    Isoflavones have a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extended SMARTS pattern for a flexible 3-aryl-1-benzopyran-4-one structure
    # The pattern allows for variability by not specifying the aryl group's exact structure
    isoflavone_pattern = Chem.MolFromSmarts('O=C1COc2cc(CC3=CC=CC=C3)ccc2C1')

    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Does not match isoflavone core structure"

    return True, "Contains isoflavone core structure (3-aryl-1-benzopyran-4-one skeleton)"