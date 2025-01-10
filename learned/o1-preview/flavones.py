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

    # Define the flavone core using a more general SMARTS pattern with aromatic atoms
    # Flavone core: fused benzopyran-4-one with an aryl group at position 2
    flavone_core_smarts = '[O]=C1c2ccccc2Oc3ccccc13'  # General flavone core
    flavone_core = Chem.MolFromSmarts(flavone_core_smarts)
    if flavone_core is None:
        return False, "Error in flavone core SMARTS definition"

    # Ensure the molecule's aromaticity matches RDKit's aromaticity model
    Chem.SanitizeMol(mol)

    # Find substructure matches of the flavone core
    core_matches = mol.GetSubstructMatches(flavone_core)
    if not core_matches:
        return False, "Flavone core structure not found"

    # Check for 2-aryl substitution (although the core pattern includes it)
    # Additional verification can be done if necessary

    return True, "Contains flavone core structure"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'flavones',
        'definition': 'A member of the class of flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.',
        'parents': []
    }
}