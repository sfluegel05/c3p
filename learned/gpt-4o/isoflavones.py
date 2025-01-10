"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Update SMARTS pattern for the 3-aryl-1-benzopyran-4-one structure
    # This pattern should capture the core structure with variations
    isoflavone_pattern = Chem.MolFromSmarts("[cC]1[cC](-c2ccccc2)[cC](=O)oc3ccccc13")

    # Check if the molecule contains the isoflavone pattern
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Does not contain the 3-aryl-chromen-4-one core structure"
    
    # It contains the essential structure of isoflavones
    return True, "Contains the 3-aryl-chromen-4-one core structure"