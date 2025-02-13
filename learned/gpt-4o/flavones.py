"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    Flavones are characterized by a 2-aryl-1-benzopyran-4-one skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a flavone core structure
    # Pattern aims to be more comprehensive and specific
    flavone_pattern = Chem.MolFromSmarts("c1c(=O)oc2cc(c(-c3ccccc3)o2)c1")

    # Check if the molecule matches the flavone pattern
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No 2-aryl-1-benzopyran-4-one skeleton found"

    # Additional verification steps can be added here

    return True, "Contains 2-aryl-1-benzopyran-4-one skeleton"