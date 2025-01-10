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
    
    # Define a generalized SMARTS pattern for the core flavanone skeleton
    # Focus on a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one structure
    flavanone_pattern = Chem.MolFromSmarts("c1ccccc1C2CC(=O)c3ccccc3O2")
    # Core feature: aromatic benzene ring, carbonyl on chromanone, 2nd aryl group attached on chromanone
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core skeleton found"
    
    # Additional checks for variability in substituents can be done if needed
    
    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"