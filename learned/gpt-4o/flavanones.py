"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    
    A flavanone is characterized by a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton
    with substitutions possible on the aromatic rings.
    
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

    # Define the flavanone core SMARTS pattern with 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one structure
    # Note: The SMARTS pattern may need refinement to better match flavanone characteristics
    flavanone_pattern = Chem.MolFromSmarts("O=C1CC[C@H](Oc2ccccc2)c2ccccc12")

    # Check for the core flavanone structure
    if not mol.HasSubstructMatch(flaonone_pattern):
        return False, "No flavanone core structure found"

    return True, "Flavanone core structure found"

# Please note: This SMARTS pattern is indicative and might need further refinement considering stereochemistry and specific substituents present in flavanones.