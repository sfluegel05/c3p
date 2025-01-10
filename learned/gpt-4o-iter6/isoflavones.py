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
    
    # Define the SMARTS pattern for 3-aryl-1-benzopyran-4-one
    # Refined pattern considering 3-aryl connectivity and possible ketone position
    isoflavone_pattern = Chem.MolFromSmarts('C1(CC(=O)Oc2ccccc2)c2ccccc2OC1')
    
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Does not match refined isoflavone core structure"

    return True, "Contains isoflavone core structure (3-aryl-1-benzopyran-4-one skeleton)"