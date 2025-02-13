"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    Flavones are characterized by a 2-aryl-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a broader and more flexible SMARTS pattern for a flavone core structure
    # This pattern includes variable attachments for characteristic flavone substitutions
    flavone_pattern = Chem.MolFromSmarts("c1cc(-c2ccoc2=O)cc1 |c:c:c:c:c:c__||c:c:o:c:3|")

    # Check if the molecule matches the flavone pattern
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No 2-aryl-1-benzopyran-4-one skeleton found"

    return True, "Contains 2-aryl-1-benzopyran-4-one skeleton"