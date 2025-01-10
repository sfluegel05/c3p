"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    Isoflavones have a 3-aryl-1-benzopyran-4-one skeleton.

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
    # Example SMARTS: O=C1C=Cc2ccccc2Oc2cccc(Cc)=c12
    isoflavone_pattern = Chem.MolFromSmarts('O=C1OCc2ccccc2C1')
    
    if mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Contains isoflavone core structure (3-aryl-1-benzopyran-4-one skeleton)"
    else:
        return False, "Does not match isoflavone core structure"