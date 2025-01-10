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
    
    # Define a SMARTS pattern for the core 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton
    # Allow substitutions on aromatic rings
    flavanone_pattern = Chem.MolFromSmarts("c1cccc2C(=O)C[C@H](Oc12)c3ccccc3")
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core skeleton found"
    
    # Further checks for optional substituents could be added here as necessary
    
    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"