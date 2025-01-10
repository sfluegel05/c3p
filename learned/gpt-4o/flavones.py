"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-aryl-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Corrected flavone core pattern, allowing for substitutions on the rings
    flavone_pattern = Chem.MolFromSmarts("c1cc(-c2cc(-c3oc(cc3=O)-c3ccccc3)c(=O)oc2)ccc1")
    
    # Check for flavone core pattern match
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains 2-aryl-1-benzopyran-4-one skeleton"
    else:
        return False, "Does not contain flavone skeleton"

# Test Cases
# This function should now be tested with known flavone and non-flavone compounds to ensure accuracy.