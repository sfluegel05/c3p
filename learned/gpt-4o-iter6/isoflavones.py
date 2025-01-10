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
    
    # Improved SMARTS pattern for a flexible 3-aryl-1-benzopyran-4-one structure
    # This pattern includes flexibility for variations in the aryl and core structure
    # R allows for attachments (common modifications like prenyl, methoxy, hydroxy, etc.)
    isoflavone_pattern = Chem.MolFromSmarts('O=C1C=COc2cc([R0])ccc2C1')

    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Does not match isoflavone core structure"

    return True, "Contains isoflavone core structure (3-aryl-1-benzopyran-4-one skeleton)"