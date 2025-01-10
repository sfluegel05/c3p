"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    
    A flavanone is characterized by a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core flavanone structure pattern
    flavanone_pattern = Chem.MolFromSmarts("O=C1CC2=C(O1)C=CC=C2")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core structure found"
    
    # Check for attached phenyl group at position 2
    phenyl_pattern = Chem.MolFromSmarts("[cR2]1[cR2][cR2][cR2][cR2][cR2]1")
    if mol.HasSubstructMatch(phenyl_pattern):
        return True, "Flavanone core structure with attached phenyl group found"
    
    return True, "Flavanone core structure found, no phenyl group attached at position 2"