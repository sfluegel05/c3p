"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns: 
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the 2-aryl-1-benzopyran-4-one core
    # The pattern allows for substitutions on the flavone rings
    flavone_pattern = Chem.MolFromSmarts("c1cc(-c2cc(-c3coc(=O)cc3)c(=O)oc2)ccc1")

    # Checking if the structure matches the flavone pattern
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains 2-aryl-1-benzopyran-4-one skeleton"
    else:
        return False, "Does not contain flavone skeleton"