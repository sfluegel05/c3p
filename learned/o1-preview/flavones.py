"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton
    and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core SMARTS pattern
    # Flavone core: 2-arylchromen-4-one skeleton
    flavone_core_smarts = 'c1cc(-c2ccccc2)c(=O)oc3ccccc13'  # Flavone core with 2-aryl substitution
    flavone_core = Chem.MolFromSmarts(flavone_core_smarts)
    if flavone_core is None:
        return False, "Error in flavone core SMARTS definition"

    # Check for substructure match of flavone core
    if not mol.HasSubstructMatch(flavone_core):
        return False, "Flavone core structure not found"

    return True, "Contains flavone core structure"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'flavones',
        'definition': 'A member of the class of flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.',
        'parents': []
    }
}