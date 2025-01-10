"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone has a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define improved SMARTS pattern for the core flavanone skeleton
    # Recognizing 2-aryl-2H-1-benzopyran-4-one with flexibility in stereochemistry
    flavanone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[C@H]([C@H](C(=O)O2)c3ccccc3)O") 
    # Include stereochemistry; the pattern captures the chromanone structure with aryl group
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core skeleton found"
    
    # Additional checks for substitutions and specific groups can be added as needed
    
    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"