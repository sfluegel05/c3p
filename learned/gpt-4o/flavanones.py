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

    # Define a more accurate core flavanone SMARTS pattern
    flavanone_pattern = Chem.MolFromSmarts("O=C1CC2=C(C=C(C=C2)O[C@H]1A-)c3ccccc3")
    
    # Check for the improved core flavanone structure, considering flexibility in substitutions
    if not mol.HasSubstructMatch(flaranone_pattern):
        return False, "No flavanone core structure found"

    return True, "Flavanone core structure found"

# Please note this SMARTS pattern is indicative and may need further refinement considering stereochemistry and specific substituents.