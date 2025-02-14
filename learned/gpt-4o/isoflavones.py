"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.

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

    # Define SMARTS pattern for isoflavone with better flexibility
    # Allowing room for substitutions, captures the core 3-aryl-1-benzopyran-4-one structure
    isoflavone_pattern = Chem.MolFromSmarts('c1cc2oc(=O)cc(c2c2ccccc2)c1')

    if mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Contains 3-aryl-1-benzopyran-4-one core structure"

    return False, "Does not contain an isoflavone core structure"