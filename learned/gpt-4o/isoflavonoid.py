"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define isoflavonoid core pattern:
    # Core benzopyran with a phenyl group attached at position 3
    isoflavonoid_pattern = Chem.MolFromSmarts("c1cc(-c2coc3c(c2)c(cc(c3)O)=O)ccc1")

    # Check if it matches the pattern
    if mol.HasSubstructMatch(isoflavonoid_pattern):
        return True, "Contains 1-benzopyran with aryl substituent at position 3, matching isoflavonoid structure"
    
    return False, "Does not match isoflavonoid structure: 1-benzopyran with aryl at position 3"